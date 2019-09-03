// GTPro Query Engine, Compressor, and Indexer.
//
// For license and copyright information, please see
// https://github.com/zjshi/gt-pro2.0/blob/master/LICENSE
//
// C++11 code formatted with
//
//     clang-format -style="{BasedOnStyle: llvm, ColumnLimit: 128}"

#if __linux__
#include <linux/version.h>
#if LINUX_VERSION_CODE > KERNEL_VERSION(2, 6, 22)
#define _MAP_POPULATE_empty
#endif
#endif

#ifdef _MAP_POPULATE_empty
#define MMAP_FLAGS (MAP_PRIVATE | MAP_POPULATE)
#else
#define MMAP_FLAGS (MAP_PRIVATE)
#endif

#include <assert.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <libgen.h>
#include <limits>
#include <mutex>
#include <regex>
#include <set>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace std;

// each 8 MB chunk will run in its own thread
// 12 MB fastq ~ 26,000 reads ~ 0.5 cpu seconds
constexpr auto segment_size = 12 * 1024 * 1024;

// allocate 15% more chunks than query threads
// that way query threads would never have to wait for I/O
// 15% is a guess
constexpr auto read_ahead_ratio = 1.15;

// The DB k-mers are 31-mers.
constexpr auto K = 31;

// 2 bits to encode each ACTG letter
constexpr auto BITS_PER_BASE = 2;

// number of bits to encode entire K-mer
constexpr auto K2 = BITS_PER_BASE * K;

constexpr uint64_t LSB = 1;
constexpr auto BASE_MASK = (LSB << BITS_PER_BASE) - LSB;
constexpr auto FULL_KMER = (LSB << K2) - LSB;

// Choose appropriately sized integer types to represent offsets into
// the database.
using LmerRange = uint64_t;
constexpr auto START_BITS = 48;
constexpr auto LEN_BITS = 64 - START_BITS;
constexpr auto MAX_START = (LSB << START_BITS) - LSB;
constexpr auto MAX_LEN = (LSB << LEN_BITS) - LSB;

// element 0, 1:  61-bp nucleotide sequence centered on SNP, for a more detailed map please see
// the comment "note on the binary representation of nucleotide sequences" below.
// element 2:  SNP coordinates consisting of species ID, major/minor allele bit, and genomic position;
// new -- the 8 MSBs are now a virtual snp id (vid) to help represent conflicting kmers

constexpr auto SNP_VID_BITS = 8;
constexpr auto SNP_MAX_VID = (LSB << SNP_VID_BITS) - LSB;
constexpr auto SNP_REAL_ID_BITS = 64 - SNP_VID_BITS;
constexpr auto SNP_MAX_REAL_ID = (LSB << SNP_REAL_ID_BITS) - LSB;

// This param is only useful for perf testing.  The setting below, not
// to exceed 64 TB of RAM, is equivalent to infinity in 2019.
constexpr auto MAX_MMAP_GB = 64 * 1024;
constexpr auto MAX_END = MAX_MMAP_GB * (LSB << 30) / 8;

size_t get_fsize(const char *filename) {
  struct stat st;
  if (stat(filename, &st) == -1) {
    // Probably file not found.
    return 0;
  }
  return st.st_size;
}

struct CodeDict {
  vector<uint8_t> code_dict;
  uint8_t *data;
  CodeDict() {
    constexpr auto CHAR_LIMIT = 1 << (sizeof(char) * 8);
    for (uint64_t c = 0; c < CHAR_LIMIT; ++c) {
      // This helps us detect non-nucleotide characters on encoding.
      code_dict.push_back(0xff);
    }
    code_dict['A'] = code_dict['a'] = 0;
    code_dict['C'] = code_dict['c'] = 1;
    code_dict['G'] = code_dict['g'] = 2;
    code_dict['T'] = code_dict['t'] = 3;
    data = code_dict.data();
  }
};
const CodeDict code_dict;

// see the comment "note on the binary representation of nucleotide sequences" below.
template <class int_type, int len> bool seq_encode(int_type *result, const char *buf) {
  result[0] = 0; // forward
  result[1] = 0; // rc
  // This loop may be unrolled by the compiler because len is a compile-time constant.
  for (int_type bitpos = 0; bitpos < BITS_PER_BASE * len; bitpos += BITS_PER_BASE) {
    const uint8_t b_code = code_dict.data[*buf++];
    if (b_code & 0xfc) {
      return true;
    }
    result[0] |= (((int_type)b_code) << bitpos);
    result[1] <<= BITS_PER_BASE;
    result[1] |= (b_code ^ 3);
  }
  return false; // success
}

uint64_t reverse_complement(uint64_t dna) {
  // The bitwise complement represents the ACTG nucleotide complement -- see seq_encode above.
  // Possibly not the fastest way to RC a kmer, but only runs during database construction.
  dna ^= FULL_KMER;
  uint64_t rc = dna & BASE_MASK;
  for (int k = 1; k < K; ++k) {
    dna >>= BITS_PER_BASE;
    rc <<= BITS_PER_BASE;
    rc |= (dna & BASE_MASK);
  }
  return rc;
}

long chrono_time() {
  using namespace chrono;
  return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

int64_t kmer_lookup_chunk(vector<uint64_t> *kmer_matches, const LmerRange *const lmer_index, const uint64_t *const mmer_bloom,
                          const uint32_t *const kmers_index, const uint64_t *const snps, const char *const window,
                          const int bytes_in_chunk, const int M2, const int M3, const char *const in_path, const long s_start) {

  const uint64_t MAX_BLOOM = (LSB << M3) - LSB;

  // Reads that contain wildcard characters ('N' or 'n') are split into
  // tokens at those wildcard characters.  Each token is processed as
  // though it were a separate read.
  constexpr int MAX_TOKEN_LENGTH = 500;
  constexpr int MIN_TOKEN_LENGTH = 31;

  // This ranges from 0 to the length of the longest read (could exceed MAX_TOKEN_LENGTH).
  int token_length = 0;
  char c = '\0';

  int n_lines = 0;

  char seq_buf[MAX_TOKEN_LENGTH];
  unordered_map<uint64_t, int> footprint;

  for (int i = 0; i < bytes_in_chunk; ++i) {

    // c is the character preceding window[i]
    if (c == '\n') {
      ++n_lines;
      // The first character on the read header must be @
      if (n_lines % 4 == 0 && window[i] != '@') {
        return -i; // here i >= 1
      }
      // Similar rule from the fastq spec.
      if (n_lines % 4 == 2 && window[i] != '+') {
        return -i;
      }
    }

    // Invariant:  The number of new line characters consumed before window[i]
    // is the value of n_lines.

    c = window[i];

    // In FASTQ format, every 4 lines define a read.  The first line is the
    // read header.  The next line is the read sequence.  We only care about
    // the read sequence, where n_lines % 4 == 1.
    if (n_lines % 4 != 1) {
      // The current line does *not* contain a read sequence.
      // Next character, please.
      continue;
    }

    // The current line contains a read sequence.  Split it into tokens at wildcard 'N'
    // characters.  Buffer current token in seq_buf[0...MAX_READ_LENGTH-1].
    const bool at_token_end = (c == '\n') || (c == 'N') || (c == 'n');
    if (!(at_token_end)) {
      // Invariant:  The current token length is token_length.
      // Only the first MAX_TOKEN_LENGTH charaters of the token are retained.
      if (token_length < MAX_TOKEN_LENGTH) {
        seq_buf[token_length] = c;
      }
      ++token_length;
      // next character, please
      continue;
    }

    // is token length within acceptable bounds?   if not, token will be dropped silently
    if (MIN_TOKEN_LENGTH <= token_length && token_length <= MAX_TOKEN_LENGTH) {

      // yes, process token
      for (int j = 0; j <= token_length - K; ++j) {

        // This kmer_tuple encodes the forward and reverse kmers at seq_buf[j...]
        uint64_t kmer_tuple[2];
        const auto malformed_dna = seq_encode<uint64_t, K>(kmer_tuple, seq_buf + j);
        if (malformed_dna) {
          return -(i - token_length);
        }
        const uint64_t kmer = min(kmer_tuple[0], kmer_tuple[1]);

        if ((mmer_bloom[(kmer & MAX_BLOOM) / 64] >> (kmer % 64)) & 1) {

          const uint32_t lmer = kmer >> M2;
          const auto range = lmer_index[lmer];
          const auto start = range >> LEN_BITS;
          const auto end = min(MAX_END, start + (range & MAX_LEN));
          for (uint64_t z = start; z < end; ++z) {
            const auto kmi = kmers_index[z];
            const auto offset = kmi & 0x1f;
            const auto snp_id = kmi >> 5;
            const auto snp_repr = snps + 3 * snp_id;
            const auto low_bits = snp_repr[0] >> (62 - (offset * BITS_PER_BASE));
            const auto high_bits = (snp_repr[1] << (offset * BITS_PER_BASE)) & FULL_KMER;
            const auto db_kmer = high_bits | low_bits;
            // The assert below is true but might be a bit slow due to reverse_complement.
            // TODO: instead of calling reverse_complement() consider using kmer_rc
            // assert(lmer == (db_kmer >> M2) || lmer == (reverse_complement(db_kmer) >> M2));
            if (kmer_tuple[0] == db_kmer || kmer_tuple[1] == db_kmer) {
              // The set of kmers that cover the SNP within the given read won't conflict
              // with each other, but may still belong to different virtual SNPs.  We need
              // to identify the real SNP they all belong to, and increment the real SNP's
              // counter just once for the read.
              const auto snp = snp_repr[2] & SNP_MAX_REAL_ID;
              if (footprint.find(snp) == footprint.end()) {
                kmer_matches->push_back(snp);
                footprint.insert({snp, 1});
              }
            }
          }
        }
      }
    }

    // clear footprint for every read instead of every token
    if (c == '\n') {
      footprint.clear();
    }

    // next token, please
    token_length = 0;
  }

  if (token_length != 0) {
    // Truncated read sequence at end of chunk.  Malformed FASTQ.
    return -bytes_in_chunk;
  }

  return (n_lines + 3) / 4;
}

char *last_read(const char *window, const int bytes_in_window) {
  // Return a pointer to the initial '@' character of the last read header
  // within the given window.  Return NULL if no such read (for example,
  // if the file does not conform to FASTQ format).

  // HOW?
  //
  // In FASTQ format, a line that starts with '@' is either a read header
  // (also called sequence identifier), or a quality string.
  //
  // Examining only the local neighborhood of an arbitrary line L that
  // starts with '@', we may quickly determine whether L contains a read header
  // or a quality string by looking at the first character of line L - 2.
  //
  // If the first character of line L - 2 is '+', then line L is a read header.
  // Otherwise, line L is a quality string and line L - 3 is a read header.
  //
  // This is not merely an heuristic.  It follows from the fastq format definition.

  // First find where in the window the last few lines begin.
  // Let line_start[3] point to the first character on the last line within the
  // window, line_start[2] same for the previous line, etc.  If fewer than 4 lines
  // (malformed fastq), let all remaining line_start's point to window[0].
  // Then compute L such that line_start[L] points to the '@' character of the
  // last read header in the window.

  // Below we use "bytes_in_window - 1" not just "bytes_in_window" because we
  // are interested in the character that comes after the newline.

  // we need bytes_in_window to be at least 2 for memrchr below
  // no read is under 4 bytes in FASTQ
  if (bytes_in_window <= 4) {
    return NULL;
  }

  char *newline = (char *)memrchr(window, '\n', bytes_in_window - 1);
  int L = -1;
  char *line_start[4];
  for (int l = 3; l >= 0; --l) {
    if (newline == NULL) {
      line_start[l] = NULL;
      continue;
    }
    line_start[l] = newline + 1;
    if (L == -1 && *(line_start[l]) == '@') {
      L = l;
    }
    newline = (char *)memrchr(window, '\n', newline - window);
  }
  if (L >= 2) {
    if (line_start[L - 2] != NULL && *(line_start[L - 2]) != '+') {
      L -= 3;
    }
  }
  if (L >= 0 && line_start[L] != NULL && *(line_start[L]) == '@') {
    return line_start[L];
  }
  return NULL;
}

bool kmer_lookup(LmerRange *lmer_index, uint64_t *mmer_bloom, uint32_t *kmers_index, uint64_t *snps, const int n_inputs,
                 char **input_paths, char *o_name, const int M2, const int M3, const int n_threads) {

  auto s_start = chrono_time();

  // In this design we have a single reader doing all the I/O (this function),
  // and up to n parallel query threads.  The reader grabs buffer space
  // for the next chunk of input, reads from the input file, launches a
  // thread to execute queries for the chunk, and repeats.  The reader may
  // block for three reasons:
  //
  //   *  I/O, waiting for bytes to be read
  //
  //   *  when the maximum number of threads are running, waiting for a
  //      thread to complete before a new thread can be launched
  //
  //   *  when all buffers are full, waiting for a thread to
  //      complete before more bytes can be read
  //
  // The reader may read-ahead some number of chunks.  That is, it will keep
  // reading even if it cannot immediately dispatch query threads, up until
  // the allocated buffer capacity is fully utilized.  This may help better
  // overlap I/O and computation in some cases.

  // Allocate space for this many input chunks.
  const int n_chunks = max(n_threads + 1, int(n_threads * read_ahead_ratio));
  vector<char> buffer(segment_size * n_chunks);
  char *buffer_addr = buffer.data();

  // This window within the buffer is ready to receive the first input chunk.
  int chunk_idx = 0;
  char *window = buffer_addr + chunk_idx * segment_size;

  // A chunk is running if and only if some query thread owns it.
  vector<bool> running;
  // Often actual_size would be smaller than segment_size, because a chunk contains
  // whole reads (any leftover bytes are carried through to the start of the next chunk).
  vector<int> actual_size;
  vector<int> chunk_channel;
  vector<int64_t> offset_in_file;
  for (int i = 0; i < n_chunks; ++i) {
    running.push_back(false);
    actual_size.push_back(0);
    offset_in_file.push_back(-1);
    chunk_channel.push_back(-1);
  }

  struct Result {
    int channel;
    const char *in_path;
    const char *o_name;
    int n_input_chunks, n_processed_chunks;
    bool finished_reading;
    bool done_with_output;
    vector<uint64_t> *p_kmer_matches;
    uint64_t chars_read;
    bool error;
    bool io_error;
    uint64_t error_pos;
    uint64_t n_reads;
    Result(int channel, const char *in_path, const char *o_name)
        : channel(channel), in_path(in_path), n_input_chunks(0), n_processed_chunks(0), finished_reading(false),
          done_with_output(false), o_name(o_name), p_kmer_matches(new vector<uint64_t>()), chars_read(0), error(false),
          error_pos(1ULL << 48), io_error(false), n_reads(0) {}
    ~Result() {
      if (p_kmer_matches) {
        delete p_kmer_matches;
        p_kmer_matches = NULL;
      }
    }
    void note_io_error() {
      error_pos = min(error_pos, chars_read);
      if (!(error)) {
        cerr << chrono_time() << ":  "
             << "[ERROR] Failed to read past position " << error_pos << " in presumed FASTQ file " << in_path << endl;
      }
      io_error = true;
      error = true;
    }
    void data_format_error(uint64_t pos = 1ULL << 48) {
      error_pos = min(error_pos, min(pos, chars_read));
      error = true;
    }
    void write_error_info() {
      assert(error);
      auto err_path = string(o_name) + "." + to_string(channel) + ".err";
      ofstream fh(err_path, ofstream::out | ofstream::binary);
      if (io_error) {
        // I/O errors are reported on stderr in realtime.  Note them in the .err file.
        fh << "[ERROR] Failed to read past position " << error_pos << " in presumed FASTQ file " << in_path << endl;
      } else {
        fh << "[ERROR] Failed to parse somewhere past position " << error_pos << " in presumed FASTQ file " << in_path << endl;
        cerr << chrono_time() << ":  "
             << "[ERROR] Failed to parse somewhere past position " << error_pos << " in presumed FASTQ file " << in_path
             << endl;
      }
      fh.close();
      delete p_kmer_matches;
      p_kmer_matches = NULL;
    }
    void write_output() {
      if (error) {
        write_error_info();
        return;
      }
      cerr << chrono_time() << ":  "
           << "[Done] searching is completed for the " << n_reads << " reads input from " << in_path << endl;
      auto out_path = string(o_name) + "." + to_string(channel) + ".tsv";
      ofstream fh(out_path, ofstream::out | ofstream::binary);
      if (p_kmer_matches->size() == 0) {
        cerr << chrono_time() << ":  "
             << "[WARNING] found zero hits for the " << n_reads << " reads input from " << in_path << endl;
      } else {
        // For each kmer output how many times it occurs in kmer_matches.
        sort(p_kmer_matches->begin(), p_kmer_matches->end());
        const uint64_t end = p_kmer_matches->size();
        uint64_t i = 0;
        uint64_t n_snps = 0;
        uint64_t n_hits = 0;
        while (i != end) {
          uint64_t j = i + 1;
          while (j != end && (*p_kmer_matches)[i] == (*p_kmer_matches)[j]) {
            ++j;
          }
          ++n_snps;
          n_hits += (j - i);
          fh << (*p_kmer_matches)[i] << '\t' << (j - i) << '\n';
          i = j;
        }
        cerr << chrono_time() << ":  "
             << "[Stats] " << n_snps << " snps, " << n_reads << " reads, " << int((((double)n_hits) / n_snps) * 100) / 100.0
             << " hits/snp, for " << in_path << endl;
      }
      fh.close();
      delete p_kmer_matches;
      p_kmer_matches = NULL;
    }
    bool pending_output() {
      if (done_with_output || !(finished_reading)) {
        return false;
      }
      {
        unique_lock<mutex> lk(mtx);
        return n_input_chunks == n_processed_chunks;
      }
    }
    void merge_kmer_matches(vector<uint64_t> &kmt, const int64_t n_reads_chunk) {
      unique_lock<mutex> lk(mtx);
      p_kmer_matches->insert(p_kmer_matches->end(), kmt.begin(), kmt.end());
      ++n_processed_chunks;
      if (n_reads_chunk >= 0) {
        n_reads += n_reads_chunk;
      }
    }

  private:
    mutex mtx;
  };

  vector<Result *> results;
  for (int i = 0; i < n_inputs; ++i) {
    results.push_back(new Result(i, input_paths[i], o_name));
  }

  int first_result_idx_not_done_with_output = 0;
  int last_result_idx_done_with_input = -1;

  uint64_t total_chars = 0;
  uint64_t total_reads = 0;

  // bytes_leftover is the number of bytes following the last complete read that is fully contained in the chunk window
  // these bytes are copied over to the next chunk window
  uint64_t bytes_leftover = 0;

  // this mutex protects "n_running_threads, "actual_size", and "running"
  mutex mtx;
  condition_variable cv;
  int n_running_threads = 0;

  int channel = 0;
  const char *in_path = input_paths[channel];
  int fd = open(in_path, O_RDONLY);

  //  Print progress update every 1 million reads.
  constexpr uint64_t PROGRESS_UPDATE_INTERVAL = 1000 * 1000;
  uint64_t total_reads_last_update = 0;

  // We may have finished_reading_all_input but some_work_remains if threads are still running
  // or if chunk data has been buffered and is waiting to be launched in new threads.
  bool some_work_remains = true;
  bool finished_reading_all_input = false;

  while (some_work_remains) {

    // the initial bytes_left_over in window represent still unprocessed input from the same file;
    // now read more bytes up to the end of the window
    ssize_t bytes_read;
    if (channel < n_inputs && results[channel]->error) {
      bytes_read = 0;
      bytes_leftover = 0;
    } else {
      bytes_read = read(fd, window + bytes_leftover, segment_size - bytes_leftover);
    }
    if (bytes_read == (ssize_t)-1) {
      assert(channel < n_inputs);
      results[channel]->note_io_error();
      continue;
    }

    const auto bytes_in_window = bytes_leftover + bytes_read;
    if (bytes_in_window == 0) {
      // we are done for this file
      close(fd);
      results[channel]->finished_reading = true;
      if (channel > last_result_idx_done_with_input) {
        last_result_idx_done_with_input = channel;
      }
      ++channel;
      if (channel == n_inputs) {
        finished_reading_all_input = true;
      } else {
        in_path = input_paths[channel];
        fd = open(in_path, O_RDONLY);
        continue;
      }
    } else {
      results[channel]->n_input_chunks++;
      chunk_channel[chunk_idx] = channel;
    }

    if (bytes_in_window > 0 && window[0] != '@') {
      results[channel]->data_format_error();
      continue;
    }

    total_chars += bytes_read;
    if (bytes_read) {
      assert(channel < n_inputs);
      results[channel]->chars_read += bytes_read;
    }

    char *last_read_addr = last_read(window, bytes_in_window);

    if (bytes_in_window > 0 && last_read_addr == NULL) {
      // TODO: Count lines here to report a more useful error message.
      assert(channel < n_inputs);
      results[channel]->data_format_error();
      continue;
    }

    const auto bytes_until_last_read = last_read_addr - window;
    const auto bytes_in_chunk = (bytes_in_window < segment_size) ? bytes_in_window : bytes_until_last_read;
    actual_size[chunk_idx] = bytes_in_chunk;
    if (bytes_in_window) {
      // the offset in file only matters if there are bytes in the chunk (otherwise the chunk won't be launched on a thread)
      assert(channel < n_inputs);
      offset_in_file[chunk_idx] = results[channel]->chars_read - bytes_in_window;
    }
    bytes_leftover = bytes_in_window - bytes_in_chunk;
    const auto leftover_addr = window + bytes_in_chunk;

    // this function will output result for an input file
    auto write_output_func = [&](const int result_idx) {
      auto &r = *results[result_idx];
      r.write_output();
      {
        unique_lock<mutex> lk(mtx); // this acquires the mutex mtx
        --n_running_threads;
        lk.unlock(); // this explicit unlock is a perf optimization of some sort, see "condition_variable" docs for C++11
        cv.notify_one();
      }
    };

    // this function will run in each query thread
    auto run_queries = [&](int my_chunk_idx) {
      int64_t n_reads;
      {
        vector<uint64_t> kmt;
        n_reads = kmer_lookup_chunk(&kmt, lmer_index, mmer_bloom, kmers_index, snps, buffer_addr + segment_size * my_chunk_idx,
                                    actual_size[my_chunk_idx], M2, M3, in_path, s_start);
        if (n_reads < 0) {
          results[chunk_channel[my_chunk_idx]]->data_format_error(offset_in_file[my_chunk_idx] - n_reads);
        }
        results[chunk_channel[my_chunk_idx]]->merge_kmer_matches(kmt, n_reads);
      }
      {
        unique_lock<mutex> lk(mtx); // this acquires the mutex mtx
        --n_running_threads;
        actual_size[my_chunk_idx] = 0;
        offset_in_file[my_chunk_idx] = -1;
        running[my_chunk_idx] = false;
        chunk_channel[my_chunk_idx] = -1;
        if (n_reads >= 0) {
          total_reads += n_reads;
        }
        lk.unlock(); // this explicit unlock is a perf optimization of some sort, see "condition_variable" docs for C++11
        cv.notify_one();
      }
    };

    // wait for thread slot or chunk slot to free up
    int ready_to_fill_chunk_idx = -1;
    while (some_work_remains && ready_to_fill_chunk_idx == -1) {
      unique_lock<mutex> lk(mtx);
      int ready_to_run_chunk_idx = -1;
      int pending_output_result_idx = -1;
      cv.wait(lk, [&] {
        // Progress update.
        if ((total_reads - total_reads_last_update) >= PROGRESS_UPDATE_INTERVAL) {
          cerr << chrono_time() << ":  " << (total_reads / 10000) / 100.0 << " million reads were scanned after "
               << (chrono_time() - s_start) / 1000 << " seconds" << endl;
          total_reads_last_update = total_reads;
        }
        // Can we emit output for some input path all of whose chunks have been processed?
        if (n_running_threads < n_threads) {
          for (int i = first_result_idx_not_done_with_output; i <= last_result_idx_done_with_input; ++i) {
            if (results[i]->pending_output()) {
              pending_output_result_idx = i;
              return true;
            }
          }
        }
        // Can we launch a thread for some chunk that's ready and waiting?
        if (n_running_threads < n_threads) {
          for (int i = 0; i < n_chunks; ++i) {
            if (!running[i] && actual_size[i] != 0) {
              ready_to_run_chunk_idx = i;
              return true;
            }
          }
        }
        // Are we totally done?  This means no more input, no more chunks to launch, and no more active threads.
        if (finished_reading_all_input && n_running_threads == 0) {
          // Here 0 == n_running_threads < n_threads, so the loop above didn't find a ready_to_run_chunk_idx.
          some_work_remains = false;
          return true;
        }
        // If we got here, we can't launch a new query thread.  Can we at least fill a new chunk window with input?
        if (!(finished_reading_all_input)) {
          for (int i = 0; i < n_chunks; ++i) {
            if (actual_size[i] == 0) {
              ready_to_fill_chunk_idx = i;
              return true;
            }
          }
        }
        // There is nothing we can do but wait (for some threads to finish).
        return false;
      });
      if (ready_to_run_chunk_idx != -1) {
        n_running_threads += 1;
        assert(actual_size[ready_to_run_chunk_idx] != 0);
        assert(!running[ready_to_run_chunk_idx]);
        running[ready_to_run_chunk_idx] = true;
        thread(run_queries, ready_to_run_chunk_idx).detach();
      } else if (pending_output_result_idx != -1) {
        n_running_threads += 1;
        auto &r = *results[pending_output_result_idx];
        r.done_with_output = true;
        thread(write_output_func, pending_output_result_idx).detach();
        while (first_result_idx_not_done_with_output < channel &&
               results[first_result_idx_not_done_with_output]->done_with_output) {
          first_result_idx_not_done_with_output += 1;
        }
      }
    }

    if (finished_reading_all_input) {
      assert(!(some_work_remains));
      assert(n_running_threads == 0);
      assert(ready_to_fill_chunk_idx == -1);
    } else {
      assert(some_work_remains);
      assert(ready_to_fill_chunk_idx != -1);
      assert(actual_size[ready_to_fill_chunk_idx] == 0);
      assert(!running[ready_to_fill_chunk_idx]);
      // move any unprocessed bytes from the old chunk window to the new chunk window
      window = buffer_addr + ready_to_fill_chunk_idx * segment_size;
      if (bytes_leftover) {
        memmove(window, leftover_addr, bytes_leftover);
      }
      chunk_idx = ready_to_fill_chunk_idx;
    }
  }
  cerr << chrono_time() << ":  " << (total_reads / 10000) / 100.0 << " million reads were scanned after "
       << (chrono_time() - s_start) / 1000 << " seconds" << endl;
  int files_with_errors = 0;
  int files_without_errors = 0;
  uint64_t reads_covered = 0;
  for (int i = 0; i < n_inputs; ++i) {
    if (results[i]->error) {
      ++files_with_errors;
    } else {
      ++files_without_errors;
      reads_covered += results[i]->n_reads;
    }
    delete results[i];
  }
  if (files_without_errors) {
    cerr << chrono_time() << ":  "
         << "Successfully processed " << files_without_errors << " input files containing " << reads_covered << " reads."
         << endl;
  }
  if (files_with_errors) {
    cerr << "Failed to process correctly " << (files_without_errors == 0 ? "all " : "") << files_with_errors << " input files."
         << endl;
  }
}

void display_usage(char *fname) {
  cout << "usage: " << fname
       << " -d <sckmerdb_path: string> -t <n_threads; int; default 1> -o <out_prefix; string; default: cur_dir/out> [-h] "
          "input1 [input2 ...]\n";
}

template <class ElementType> struct DBIndex {

  string filename;
  bool loaded_or_mmapped;

  DBIndex(const string &filename, const uint64_t expected_element_count = 0)
      : filename(filename), mmapped_data(NULL), loaded_or_mmapped(false), expected_element_count(expected_element_count),
        fd(-1), filesize(0) {}

  ElementType *address() {
    if (mmapped_data) {
      return mmapped_data;
    }
    assert(elements.size() > 0);
    return &(elements[0]);
  }

  uint64_t elementCount() {
    if (mmapped_data) {
      return filesize / sizeof(ElementType);
    } else {
      return elements.size();
    }
  }

  vector<ElementType> *getElementsVector() { return &(elements); }

  // If file exists and nonempty, preload or mmap depending on argument, and return false.
  // If file is missing or empty, allocate space in elements array and return true.
  bool mmap_or_load(const bool preload = false) {
    assert(!(loaded_or_mmapped));
    filesize = get_fsize(filename.c_str());
    if (filesize) {
      if (!(expected_element_count)) {
        expected_element_count = filesize / sizeof(ElementType);
      }
      assert(filesize == expected_element_count * sizeof(ElementType));
      if (preload) {
        load();
      } else {
        MMAP();
      }
    }
    if (loaded_or_mmapped) {
      // Does not need to be recomputed.
      return false;
    }
    cerr << chrono_time() << ":  Failed to MMAP or preload " << filename
         << ".  This is fine, but init will be slower as we recreate this file." << endl;
    elements.resize(expected_element_count);
    // needs to be recomputed
    return true;
  }

  void save() {
    assert(!(loaded_or_mmapped));
    auto l_start = chrono_time();
    FILE *dbout = fopen(filename.c_str(), "wb");
    assert(dbout);
    const auto saved_element_count = fwrite(address(), sizeof(ElementType), elements.size(), dbout);
    fclose(dbout);
    assert(saved_element_count == elements.size());
    cerr << chrono_time() << ":  Done writing " << filename << ". That took " << (chrono_time() - l_start) / 1000
         << " more seconds." << endl;
  }

  ~DBIndex() {
    if (fd != -1) {
      assert(mmapped_data);
      assert(loaded_or_mmapped);
      int rc = munmap(mmapped_data, filesize);
      assert(rc == 0);
      close(fd);
    }
  }

private:
  vector<ElementType> elements;
  ElementType *mmapped_data;
  uint64_t expected_element_count;
  int fd;
  uint64_t filesize;

  void load() {
    FILE *dbin = fopen(filename.c_str(), "rb");
    if (dbin) {
      cerr << chrono_time() << ":  Loading " << filename << endl;
      elements.resize(expected_element_count);
      const auto loaded_element_count = fread(address(), sizeof(ElementType), expected_element_count, dbin);
      if (loaded_element_count == expected_element_count) {
        cerr << chrono_time() << ":  Loaded " << filename << endl;
        loaded_or_mmapped = true;
      }
      fclose(dbin);
    }
  }

  void MMAP() {
    fd = open(filename.c_str(), O_RDONLY, 0);
    if (fd != -1) {
      cerr << chrono_time() << ":  MMAPPING " << filename << endl;
      auto mmappedData = (ElementType *)mmap(NULL, filesize, PROT_READ, MMAP_FLAGS, fd, 0);
      if (mmappedData != MAP_FAILED) {
        mmapped_data = mmappedData;
        loaded_or_mmapped = true;
        // cerr << chrono_time() << ":  MMAPPED " << filename << endl;
      }
    }
  }
};

struct SNPSeq {
  uint64_t low_64;
  uint64_t high_64;
};

int main(int argc, char **argv) {

  extern char *optarg;
  extern int optind;

  bool dbflag = false;
  bool inflag = false;

  char *fname = argv[0];
  char *db_path = (char *)"";
  char *oname = (char *)"./out";

  // Number of bits in the prefix part of the K-mer (also called L-mer,
  // even though it might not correspond to an exact number of bases).
  // Override with command line -l parameter.
  //
  // This has a substantial effect on memory use.  Rule of thumb for
  // perf is L2 >= K2 - M3.  However, that rule may be broken in order
  // to reduce RAM use and eliminate I/O which is even worse for perf.
  auto L2 = 29;

  // Number of bits in the MMER_BLOOM index.  This has a substantial effect
  // on memory use.  Rule of thumb for perf is M3 >= 4 + log2(DB cardinality).
  // Override with command line -m parameter.
  auto M3 = 36;

  int n_threads = 1;

  auto preload = false;

  int opt;
  while ((opt = getopt(argc, argv, "pl:m:d:t:o:h")) != -1) {
    switch (opt) {
    case 'd':
      dbflag = true;
      db_path = optarg;
      break;
    case 't':
      n_threads = stoi(optarg);
      break;
    case 'o':
      oname = optarg;
      break;
    case 'l':
      L2 = stoi(optarg);
      break;
    case 'm':
      M3 = stoi(optarg);
      break;
    case 'p':
      preload = true;
      break;
    case 'h':
    case '?':
      display_usage(fname);
      exit(1);
    }
  }

  const auto M2 = K2 - L2;

  assert(L2 > 0);
  assert(L2 <= 32);
  assert(M2 > 0);
  assert(M2 < 64);
  assert(M3 > 0);
  assert(M3 < 64);

  const auto LMER_MASK = (LSB << L2) - LSB;
  const auto MMER_MASK = (LSB << M2) - LSB;
  const auto MAX_BLOOM = (LSB << M3) - LSB;

  cout << fname << '\t' << db_path << '\t' << n_threads << "\t" << (preload ? "preload" : "mmap") << "\t" << L2 << "\t" << M3
       << endl;

  if (!dbflag) {
    cout << "missing argument: -d <sckmerdb_path: string>\n";
    display_usage(fname);
    exit(1);
  }

  int in_pos = optind;
  if (optind == argc) {
    cout << "missing argument: input (>1)\n";
    display_usage(fname);
    exit(1);
  }

  auto l_start = chrono_time();
  cerr << chrono_time() << ":  "
       << "Starting to load DB: " << db_path << endl;

  uint64_t db_filesize = get_fsize(db_path);

  string dbbase = string(basename(db_path));
  dbbase = regex_replace(dbbase, regex("\\.bin$"), "");
  dbbase = regex_replace(dbbase, regex("\\."), "_");

  if (preload) { // force preload, probably unnecessary and wasteful
    cerr << chrono_time() << ":  DB indexes will be preloaded." << endl;
  }

  // The input (un-optimized) DB is a sequence of 56-bit snp followed by 1-bit forward/rc indicator,
  // then 7 bit offset of SNP within kmer, then 64-bit kmer.  The 56-bit snp encodes the species id,
  // major/minor allele, and genomic position.  From that we build the "optimized" DBs
  // below.  The first one, db_snps, lists the unique SNPs in arbitrary order; and
  // for each SNP in addition to the 56-bits mentioned above it also shows the
  // sequence of 61bp centered on the SNP inferred from all kmers in the original DB.
  DBIndex<uint64_t> db_snps(dbbase + "_optimized_db_snps.bin");
  const bool recompute_snps = db_snps.mmap_or_load(preload);

  // This encodes the list of all kmers, sorted in increasing order.  Each kmer is represented
  // not by the 62 bits of its 31-bp nucleotide sequence but rather by 27-bits that represent
  // an index into the db_snps table above, and 5 bits representing the SNP position within
  // the kmer;  possibly using the kmer's reverse complement instead of the kmer.
  DBIndex<uint32_t> db_kmer_index(dbbase + "_optimized_db_kmer_index.bin");
  const bool recompute_kmer_index = db_kmer_index.mmap_or_load(preload);

  // Bit vector with one presence/absence bit for every possible M3-bit kmer suffix (the M3
  // LSBs of a kmer's nucleotide sequence).
  DBIndex<uint64_t> db_mmer_bloom(dbbase + "_optimized_db_mmer_bloom_" + to_string(M3) + ".bin", (1 + MAX_BLOOM) / 64);
  const bool recompute_mmer_bloom = db_mmer_bloom.mmap_or_load(preload);

  // For every kmer in the original DB, the most-signifficant L2 bits of the kmer's nucleotide sequence
  // are called that kmer's lmer.  Kmers that share the same lmer occupy a range of consecutive
  // positions in the kmer_index, and that range is lmer_index[lmer].
  DBIndex<LmerRange> db_lmer_index(dbbase + "_optimized_db_lmer_index_" + to_string(L2) + ".bin", 1 + LMER_MASK);
  const bool recompute_lmer_index = db_lmer_index.mmap_or_load(preload);
  LmerRange *lmer_index = db_lmer_index.address();

  assert(recompute_kmer_index == recompute_snps &&
         "Please delete all of the optimized DB bin files before recomputing any of them.");

  const bool recompute_everything = recompute_kmer_index || recompute_snps || recompute_lmer_index || recompute_mmer_bloom;

  int fd = -1;
  uint64_t *db_data = NULL;

  auto t_last_progress_update = chrono_time();

  if (recompute_kmer_index || recompute_snps) {
    const auto MAX_SNPS = LSB << 27;               // sorry this is a hard constraint at the moment
    const auto MAX_KMERS = LSB << 40;              // ~1 trillon, up to 5.7 billion tested
    const uint64_t MAX_INPUT_DB_KMERS = MAX_KMERS; // this is a test/debug/developer-only feature
    fd = open(db_path, O_RDONLY, 0);
    assert(fd != -1);
    db_data = (uint64_t *)mmap(NULL, db_filesize, PROT_READ, MMAP_FLAGS, fd, 0);
    assert(db_data != MAP_FAILED);
    unordered_map<uint64_t, uint32_t> snps_map;
    set<uint64_t> virtual_snps;
    vector<tuple<uint64_t, uint64_t>> snps_known_bits; // not persisted, just for integrity checking during construction
    auto &kmer_index = *db_kmer_index.getElementsVector();
    auto &snps = *db_snps.getElementsVector();
    for (uint64_t end = 0; (end < min(MAX_INPUT_DB_KMERS, db_filesize / 8)) && (kmer_index.size() < MAX_KMERS); end += 2) {
      if (((end + 2) % (20 * 1000 * 1000)) == 0 && (chrono_time() - t_last_progress_update >= 10 * 1000)) {
        // print progress update every 10 million kmers / but not more often than every 10 seconds
        t_last_progress_update = chrono_time();
        const auto t_elapsed = t_last_progress_update - l_start;
        const auto fraction_complete = double(end + 2) / (db_filesize / 8);
        const auto t_remaining = (t_elapsed / fraction_complete) * (1.0 - fraction_complete);
        cerr << t_last_progress_update << ":  Loaded almost " << ((end + 2) / (2000 * 1000)) << " million input DB kmers ("
             << int(fraction_complete * 1000) / 10.0 << " percent).  Expect to complete processing in "
             << int(t_remaining / (60 * 100)) / 10.0 << " more minutes." << endl;
      }
      const auto snp_with_offset = db_data[end];
      const auto rc = snp_with_offset & 0x80; // 0 -> forward, 0x80 -> reverse complement
      const auto orig_snp = snp_with_offset >> 8;
      auto offset = snp_with_offset & 0x7f;
      const auto db_kmer = db_data[end + 1];
      const auto db_kmer_rc = reverse_complement(db_kmer);
      if (db_kmer > db_kmer_rc) {
        // We have already seen and encoded db_kmer_rc.
        continue;
      }
      // Encode the kmer that is flagged as "forward" in the input DB, to represent the pair.
      const auto kmer = rc ? db_kmer_rc : db_kmer;
      offset = rc ? 30 - offset : offset;
      assert(orig_snp <= SNP_MAX_REAL_ID); // this should leave 8 bits at the top for virtual SNP id
                                           // 96-97% of the time, vid=0 works
                                           // the rest of the time, we create virtual SNP ids to represent all conflicting kmers
      for (uint64_t vid = 0; vid <= SNP_MAX_VID; ++vid) {
        const auto snp = orig_snp | (vid << SNP_REAL_ID_BITS);
        const auto map_entry = snps_map.find(snp);
        uint32_t snp_id;
        if (map_entry == snps_map.end()) {
          snp_id = snps_map.size();
          snps_map.insert({snp, snp_id});
          if (snp_id == MAX_SNPS) {
            cerr << chrono_time() << ":  Too many SNPs in database.  Only up to " << MAX_SNPS
                 << " are supported, and it's best to leave 10% margin for 'virtual' SNPs." << endl;
          }
          if (snp_id >= MAX_SNPS) {
            break; // break out of vid loop
          }
          snps.push_back(0);
          snps.push_back(0);
          snps.push_back(snp);
          snps_known_bits.push_back(tuple<uint64_t, uint64_t>(0, 0));
        } else {
          snp_id = map_entry->second;
        }
        assert(0 <= offset && offset < K && offset <= 31);
        auto *snp_repr = &(snps[3 * snp_id]);
        auto &snp_known_bits_mask = snps_known_bits[snp_id];
        //
        // The SNP position "offset" divides the kmer binary representation into
        // low_bits and high_bits, as follows:
        //
        //     kmer nucleotide 0 --> kmer binary bits 0, 1                                  |
        //     kmer nucleotide 1 --> kmer binary bits 2, 3                                  |  "low bits"
        //     ...                                                                          |
        //     kmer nucleotide offset --> kmer binary bits 2 * offset, 2 * offset + 1      >*<  SNP position
        //     ...                                                                          |
        //     kmer nucleotide 29 --> kmer binary bits 59, 60                               |  "high bits"
        //     kmer nucleotide 30 --> kmer binary bits 61, 62                               |
        //
        // That is, the kmer nucleotides at positions offset, offset + 1, ... are
        // encoded by the kmer binary bits at 2 * offset, 2 * offset + 1, ..., 62.
        // These so called "high_bits" match the LSBs of snp_repr[1].
        //
        // Similarly, the "low_bits" of the kmer binary representation, bits
        // 0, 1, ..., 2 * offset, 2 * offset + 1, form the MSBs of snp_repr[0].
        //
        // Obligatory ASCII art diagram below.
        /*

                                  <---------- memory addreses increase right to left ------------+

                                            62-2*offset
                                               kmer
                                             high bits       2*offset+2
                                         +<-------------->+     kmer
                                         |                |   low_bits
                                         |             +<------------------>+
                                         |             |  |                 |
                                      MSB|             |  |                 |
                                      +-------------------------------------+
                                      |00|abcd  ... ijk|XY|uvw  ...      rst|   kmer
                 MSB                  +-------------------------------------+
                 +----------------------------------------+                 |
        snp[1]   |00|???  ...          ??|abcd  ... ijk|XY|                 |
                 +----------------------------------------+                 |                 LSB
                                                       +----------------------------------------+
                                                       |XY|uvw  ...      rst|????  ...     ??|00|   snp[0]
                                                       +----------------------------------------+
                                                       |  |                 |  62-2*offset bits |
                                                       +--+                 <------------------->
                                                     SNP bits
                                                       "XY"



                  snp[0]:            64 bits, 2 LSBs unused (00)

                  snp[1], kmer:      64 bits, 2 MSBs unused (00)

                  allele XY:         MSBs of snp[0] =

                                     LSBs of snp[1] =

                                     bits offset*2, offset*2+1 of kmer



        */
        // Note the nucleotide at the SNP position is redundantly represented,
        // being included in both the high_bits and the low_bits.  Thus,
        // the 2 LSBs of snp_repr[1] always equal the 2 MSBs of snp_repr[0].
        // This is intended as a correctnes check.
        //
        // The 2 LSBs of snp_repr[0] and the 2 MSBs of snp_repr[1] are unused,
        // and empty for future extensions.
        //
        // As we construct snp_repr from kmers, we note which of the snp_repr bits
        // have been initialized so far, because future kmers for the snp must
        // agree on those bits with past kmers.  This too is a correctness check.
        // We do not persist these coverage masks to file.
        //
        // Finally, after we have constructed the optimized DB, we use it to reconstruct
        // from that the original DB, and compare to the original.  That is the final
        // and most definitive correctness check.
        //
        static_assert(sizeof(kmer) * 8 == 64, "Use uint64_t for kmers, please.");
        const auto low_bits = kmer << (62 - (offset * BITS_PER_BASE));
        const auto high_bits = kmer >> (offset * BITS_PER_BASE);
        assert(((low_bits >> 62) == (high_bits & 0x3)) && "SNP position differs in two supposedly redundant representations.");
        auto &snp_mask_0 = get<0>(snp_known_bits_mask);
        auto &snp_mask_1 = get<1>(snp_known_bits_mask);
        // We've added information to the snp_repr.  Extend the coverage masks.
        const auto kmer_mask_0 = (FULL_KMER << (62 - (offset * BITS_PER_BASE)));
        const auto kmer_mask_1 = (FULL_KMER >> (offset * BITS_PER_BASE));
        const auto mask_0 = snp_mask_0 & kmer_mask_0;
        const auto mask_1 = snp_mask_1 & kmer_mask_1;
        if (((mask_0 & snp_repr[0]) != (mask_0 & low_bits)) || ((mask_1 & snp_repr[1]) != (mask_1 & high_bits))) {
          virtual_snps.insert(snp);
          if (vid == SNP_MAX_VID) {
            // This really shouldn't happen.
            cerr << "ERROR:  SNP " << snp << " covered by " << (SNP_MAX_VID + 1) << " conflicting kmers." << endl;
            assert(false);
          }
          // continue to attempt next vid
        } else {
          snp_repr[0] |= low_bits;
          snp_repr[1] |= high_bits;
          // We've added information to the snp_repr.  Extend the coverage masks.
          snp_mask_0 |= kmer_mask_0;
          snp_mask_1 |= kmer_mask_1;
          const auto kmer_repr = (snp_id << 5) | offset;
          kmer_index.push_back(kmer_repr);
          break; // break out of VID loop as we've successfully encoded the kmer into this virtual SNP id
        }
      }
    }
    if (snps_map.size() >= MAX_SNPS) {
      cerr << chrono_time() << ":  ERROR:  SNP count exceeds maximum by " << (snps_map.size() - MAX_SNPS) << " SNPs ("
           << int((snps_map.size() - MAX_SNPS) * 1000.0 / snps_map.size()) / 10.0
           << " percent).  A margin of 5-10% is recommended." << endl;
      // We cannot ensure DB integrity with SNP overflow.  It's too complicated with virtual SNPs.
      assert(false && "too many SNPs");
    }
    cerr << chrono_time() << ":  Compressed DB contains " << (snps.size() / 3) << " snps." << endl;
    cerr << chrono_time() << ":  There were " << virtual_snps.size() << " virtual snps." << endl;
    cerr << chrono_time() << ":  Validating optimized DB against original DB." << endl;
    assert((snps.size() / 3) == snps_map.size());
    for (uint64_t end = 0, last_kmi = 0; (end < min(MAX_INPUT_DB_KMERS, db_filesize / 8)) && (last_kmi < MAX_KMERS); end += 2) {
      const auto db_snp_with_offset = db_data[end];
      const auto db_snp = db_snp_with_offset >> 8;
      const auto db_offset = db_snp_with_offset & 0x1f;
      const auto rc = db_snp_with_offset & 0x80; // 0 -> forward, 0x80 -> reverse complement
      const auto map_entry = snps_map.find(db_snp);
      assert(map_entry != snps_map.end());
      const auto db_virtual_snp_id = map_entry->second;
      assert(db_virtual_snp_id < MAX_SNPS);
      const auto db_kmer = db_data[end + 1];
      const auto db_kmer_rc = reverse_complement(db_kmer);
      if (db_kmer_rc < db_kmer) {
        // We have already decoded db_kmer_rc.
        continue;
      }
      const auto kmi = kmer_index[last_kmi++];
      const auto offset = kmi & 0x1f;
      const auto virtual_snp_id = kmi >> 5;
      assert(0 <= offset && offset <= 31 && offset < K);
      assert(0 <= virtual_snp_id && virtual_snp_id <= (snps.size() / 3));
      const auto snp = snps[3 * virtual_snp_id + 2] & SNP_MAX_REAL_ID;
      assert(snp == db_snp);
      assert(offset == (rc ? 30 - db_offset : db_offset));
      const auto *snp_repr = &(snps[3 * virtual_snp_id]);
      const auto low_bits = snp_repr[0] >> (62 - (offset * BITS_PER_BASE));
      const auto high_bits = (snp_repr[1] << (offset * BITS_PER_BASE)) & FULL_KMER;
      assert(((snp_repr[0] >> 62) == (snp_repr[1] & 0x3)) &&
             "SNP position differs in two supposedly redundant representations.");
      const auto kmer = high_bits | low_bits;
      if (kmer != db_kmer && kmer != db_kmer_rc) {
        cerr << chrono_time() << ":  ERROR:  Mismatch between original and reconstructed kmer at input DB position " << end
             << endl;
        cerr << chrono_time() << ":  ERROR:       original " << bitset<64>(db_kmer) << endl;
        cerr << chrono_time() << ":  ERROR:    original_rc " << bitset<64>(db_kmer_rc) << endl;
        cerr << chrono_time() << ":  ERROR:  reconstructed " << bitset<64>(kmer) << endl;
        cerr << chrono_time() << ":  ERROR:  snp " << db_snp << " vid "
             << ((snps[3 * virtual_snp_id + 2] ^ snp) >> SNP_REAL_ID_BITS) << " offset " << offset << endl;
        assert(false);
      }
    }
  }

  if (recompute_lmer_index || recompute_mmer_bloom) {
    cerr << chrono_time() << ":  Recomputing bloom index and/or filter." << endl;
    uint64_t start = 0;
    uint64_t last_lmer;
    const auto kmer_index = db_kmer_index.address();
    const auto kmer_count = db_kmer_index.elementCount();
    const auto snps = db_snps.address();
    const auto snps_count = db_snps.elementCount() / 3;
    auto mmer_bloom = db_mmer_bloom.address();
    t_last_progress_update = chrono_time();
    for (uint64_t end = 0; end < kmer_count; ++end) {
      if (((end + 1) % (10 * 1000 * 1000)) == 0 && (chrono_time() - t_last_progress_update >= 10 * 1000)) {
        // print progress update every 10 million kmers / but not more often than every 10 seconds
        t_last_progress_update = chrono_time();
        const auto t_elapsed = t_last_progress_update - l_start;
        const auto fraction_complete = double(end + 1) / kmer_count;
        const auto t_remaining = (t_elapsed / fraction_complete) * (1.0 - fraction_complete);
        cerr << t_last_progress_update << ":  Processed almost " << ((end + 1) / (1000 * 1000)) << " million kmers ("
             << int(fraction_complete * 1000) / 10.0 << " percent).  Expect to complete processing in "
             << int(t_remaining / (60 * 100)) / 10.0 << " more minutes." << endl;
      }
      const auto kmi = kmer_index[end];
      const auto offset = kmi & 0x1f;
      const auto snp_id = kmi >> 5;
      assert(0 <= offset && offset <= 31 && offset < K);
      assert(0 <= snp_id && snp_id <= snps_count);
      const auto *snp_repr = &(snps[3 * snp_id]);
      const auto low_bits = snp_repr[0] >> (62 - (offset * BITS_PER_BASE));
      const auto high_bits = (snp_repr[1] << (offset * BITS_PER_BASE)) & FULL_KMER;
      assert(((snp_repr[0] >> 62) == (snp_repr[1] & 0x3)) &&
             "SNP position differs in two supposedly redundant representations.");
      auto kmer = high_bits | low_bits;
      const auto kmer_rc = reverse_complement(kmer);
      if (kmer_rc < kmer) {
        // remember only the smaller of the kmer pair is in the DB, index, and bloom filter
        kmer = kmer_rc;
      }
      const auto lmer = kmer >> M2;
      if (recompute_lmer_index) {
        if ((end > 0) && (lmer != last_lmer)) {
          start = end;
        }
        // Invariant:  The data loaded so far for lmer reside at kmer_index[start...]
        assert(start <= MAX_START);
        const auto len = end - start + 1;
        assert(len < MAX_LEN);
        assert(lmer <= LMER_MASK);
        lmer_index[lmer] = (start << LEN_BITS) | len;
        last_lmer = lmer;
      }
      if (recompute_mmer_bloom) {
        const uint64_t bloom_index = kmer & MAX_BLOOM;
        mmer_bloom[bloom_index / 64] |= ((uint64_t)1) << (bloom_index % 64);
      }
    }
  }

  if (recompute_kmer_index) {
    db_snps.save();
    db_kmer_index.save();
  }

  if (recompute_mmer_bloom) {
    db_mmer_bloom.save();
  }

  if (recompute_lmer_index) {
    db_lmer_index.save();
  }

  cerr << chrono_time() << ":  Done with init for optimized DB with " << db_kmer_index.elementCount() << " kmers.  That took "
       << (chrono_time() - l_start) / 1000 << " seconds." << endl;

  l_start = chrono_time();

  kmer_lookup(lmer_index, db_mmer_bloom.address(), db_kmer_index.address(), db_snps.address(), argc - optind, argv + optind,
              oname, M2, M3, n_threads);

  if (fd != -1 && db_data != NULL) {
    int rc = munmap(db_data, db_filesize);
    assert(rc == 0);
    close(fd);
  }

  cerr << chrono_time() << ":  "
       << " Totally done: " << (chrono_time() - l_start) / 1000 << " seconds elapsed processing reads, after DB was loaded."
       << endl;

  return 0;
}
