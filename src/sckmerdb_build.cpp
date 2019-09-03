// SNP-covering K-mer database binary encoder.
//
// For license and copyright information, please see
// https://github.com/zjshi/gt-pro2.0/blob/master/LICENSE
//
// C++11 code formatted with
//
//     clang-format -style="{BasedOnStyle: llvm, ColumnLimit: 128}"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mutex>
#include <regex>
#include <thread>
#include <unordered_map>
#include <vector>

#include <assert.h>
#include <fcntl.h>
#include <libgen.h>
#include <unistd.h>

using namespace std;

// global variable declaration starts here
constexpr auto K = 31;

// Data layout:
//
// On little-endian achitectures, the second field of a tuple is stored at a lower memory address
// than the first field.  This layout is the exact opposite of what you get from a struct.
//
// The rationale behind the tuple layout is this:  The lexicographic sort order on tuples is the
// same as the numeric sort order when the tuple's bits are interpreted as one giant integer.
//
// In our application, the tuples represent a kmer and its snp.  We like to sort by the kmer first,
// so the kmer would be the first element of the tuple, will take the higher addresses in
// memory, and will come *after* the SNP when tuples are persisted to a file.
//
using KmerData = tuple<uint64_t, uint64_t>;

// output file path
constexpr auto OUT_PATH = "/dev/stdout";

// get time elapsed since when it all began in milliseconds.
long chrono_time() {
  using namespace chrono;
  return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

struct CodeDict {
  vector<int8_t> code_dict;
  int8_t *data;
  CodeDict() {
    constexpr auto CHAR_LIMIT = 1 << (sizeof(char) * 8);
    for (uint64_t c = 0; c < CHAR_LIMIT; ++c) {
      // This helps us detect non-nucleotide characters on encoding.
      code_dict.push_back(-1);
    }
    code_dict['A'] = code_dict['a'] = 0;
    code_dict['C'] = code_dict['c'] = 1;
    code_dict['G'] = code_dict['g'] = 2;
    code_dict['T'] = code_dict['t'] = 3;
    data = code_dict.data();
  }
};

const CodeDict code_dict;

template <int len> uint64_t seq_encode(const char *buf) {
  constexpr int BITS_PER_BASE = 2;
  uint64_t seq_code = 0;
  for (int i = 0; i < len; ++i) {
    const auto b_code = code_dict.data[buf[i]];
    assert(b_code != -1 && "Trying to encode a character that is not ACTG or actg.");
    seq_code |= (((uint64_t)b_code) << (BITS_PER_BASE * i));
  }
  return seq_code;
}

string profile_path(string fname) {
  // fname must be DDDDDD.sckmer_allowed.tsv
  // result would be DDDDDD.sckmer_profiles.tsv
  string base = basename(const_cast<char *>(fname.c_str()));
  string dir = dirname(const_cast<char *>(fname.c_str()));
  string altbase =
      regex_replace(base, regex("\\b([1-9][0-9][0-9][0-9][0-9][0-9]).sckmer_allowed.tsv$"), "$1.sckmer_profiles.tsv");
  if (altbase.length() != strlen("DDDDDD.sckmer_profiles.tsv") || base == altbase) {
    cerr << "Malformed argument: " << base << endl;
    cerr << "Required format: DDDDDD.sckmer_allowed.tsv" << endl;
    cerr << "Additional details: " << altbase << endl;
    assert(false);
  }
  if (dir == ".") {
    return altbase;
  }
  return dir + "/" + altbase;
}

inline uint64_t get_kmer(const KmerData &kmer_data) { return get<0>(kmer_data); }

inline uint64_t get_snp_with_rc_and_offset(const KmerData &kmer_data) { return get<1>(kmer_data); }

inline uint64_t get_snp(const KmerData &kmer_data) { return get<1>(kmer_data) >> 8; }

// INPUT:  A sckmers.tsv file with content like:
//
// ...
//
// 343 25  TGGATACCACGGCGCAAGAGCACGCACAGAA TGGATACCACGGCGCAAGAGCACGCGCAGAA \
//   TTCTGTGCGTGCTCTTGCGCCGTGGTATCCA TTCTGCGCGTGCTCTTGCGCCGTGGTATCCA 1   11  259564  8   3
//
// 343 24  GGATACCACGGCGCAAGAGCACGCACAGAAT GGATACCACGGCGCAAGAGCACGCGCAGAAT \
//   ATTCTGTGCGTGCTCTTGCGCCGTGGTATCC ATTCTGCGCGTGCTCTTGCGCCGTGGTATCC 1   11  259564  8   3
//
// ...
//
// where each line contains, in this order:
//
//      $1 - genomic position of SNP
//      $2 - zero-based offset of SNP in forward 31-mer
//      $3 - forward 31-mer covering the SNP's major allele
//      $4 - forward 31-mer covering the SNP's minor allele
//      $5 - reverse complement of $3
//      $6 - reverse complement of $4
//      $7 - unused
//      $8 - unused
//      $9 - unused
//      $10 - 6 decimal digit species ID
//      $11 - unused
//      $12 - unused
//
//  a "snp coordinate" is obtained by concatenating the following character strings and
//  converting the result back to an integer:
//
//      * the 6 decimal digit species id (input column $10),
//
//      * a single digit allele type:
//           (0 = major allele, matching input $3 and $5,
//            1 = minor allele, matching input $4 and $6),
//
//      * genomic position of SNP (input column $1)
//

void bit_load(string f_allowed, vector<KmerData> &result, vector<KmerData> &profiles) {

  auto t_start = chrono_time();

  string f_profile = profile_path(f_allowed);

  // String addition so hopefully this << is atoimic.
  cerr << (string("Loading files ") + f_allowed + " and " + f_profile + "\n");

  FILE *fp = fopen(f_profile.c_str(), "r");
  if (fp == NULL) {
    cerr << "Trouble opening file " << f_profile << "." << endl;
    assert(false);
  }

  FILE *fa = fopen(f_allowed.c_str(), "r");
  if (fa == NULL) {
    cerr << "Trouble opening file " << f_allowed << "." << endl;
    assert(false);
  }

  struct Column {
    char *c_str;
    size_t size;
    ssize_t len;
    int sep;
    const char *name;
    Column(const char *name, int separator) : c_str(NULL), size(0), len(0), sep(separator), name(name){};
    ~Column() {
      if (c_str) {
        free(c_str);
      }
    }
    inline void read_from_file(FILE *fp) {
      // getdelim will (re)allocate memory for c_str as needed.
      len = getdelim(&c_str, &size, sep, fp);
      // squash delimiter if present
      if (len > 0 && c_str[len - 1] == sep) {
        c_str[--len] = '\0';
      }
    }
    inline bool is_numeric() {
      char *s = c_str;
      while (isdigit(*s)) {
        ++s;
      }
      return *s == '\0';
    }
    inline bool does_not_contain_wildcard() {
      char *s = c_str;
      while (*s && *s != 'n' && *s != 'N') {
        ++s;
      }
      return *s == '\0';
    }
    inline bool is_snp_id() { return ((len >= 8) && (len <= 16) && is_numeric() && ((c_str[6] == '0') || (c_str[6] == '1'))); }
  };

  // columns in the DDDDDD.sckmer_allowed.tsv file.
  Column allowed_cols[] = {Column("kmer", '\t'), Column("snp_id", '\n')};
  auto &ac_kmer = allowed_cols[0];
  auto &ac_snp = allowed_cols[1];

  uint64_t allowed_count = 0;
  while (true) {
    bool no_column_is_empty = true;
    for (auto &c : allowed_cols) {
      c.read_from_file(fa);
      no_column_is_empty = no_column_is_empty && c.len > 0;
    }
    if (allowed_cols[0].len < 0) {
      // end of file
      break;
    }
    assert(no_column_is_empty && "Empty columnns are not allowed.");
    assert(ac_kmer.len == K && "Kmer column must contain exactly K characters");
    assert(ac_snp.is_snp_id() &&
           "SNP id column needs to be decimal with at 8...16 digits and the character at position 6 must be '0' or '1'.");
    auto kmer = seq_encode<K>(ac_kmer.c_str);
    auto snp = strtoull(ac_snp.c_str, NULL, 10);
    result.push_back(KmerData(kmer, snp)); // the snp will be amended later with more information
    ++allowed_count;
  }
  fclose(fa);
  auto end = result.end();
  auto start = end - allowed_count;
  sort(start, end);

  // columns in the DDDDDD.sckmer_profiles.tsv file.
  Column profile_cols[] = {Column("snp_genomic_position", '\t'),
                           Column("snp_offset_in_forward_kmer", '\t'),
                           Column("forward_kmer_major_allele", '\t'),
                           Column("forward_kmer_minor_allele", '\t'),
                           Column("reverse_complement_kmer_major_allele", '\t'),
                           Column("reverse_complement_kmer_minor_allele", '\t'),
                           Column("unused", '\t'),
                           Column("unused", '\t'),
                           Column("unused", '\t'),
                           Column("genome_id", '\t'),
                           Column("unused", '\t'),
                           Column("unused", '\n')};

  auto kmer_columns = {2, 3, 4, 5};

  auto &c_snp_genomic_position = profile_cols[0];
  auto &c_snp_offset_in_forward_kmer = profile_cols[1];
  auto &c_genome_id = profile_cols[9];

  profiles.clear();
  while (true) {
    bool no_column_is_empty = true;
    for (auto &c : profile_cols) {
      c.read_from_file(fp);
      no_column_is_empty = no_column_is_empty && c.len > 0;
    }
    if (profile_cols[0].len < 0) {
      // end of file
      break;
    }
    assert(no_column_is_empty && "Empty columnns are not allowed.");
    for (auto kc : kmer_columns) {
      assert(profile_cols[kc].len == K && "Kmer column must contain exactly K characters");
    }
    assert(c_snp_genomic_position.len <= 9 && c_snp_genomic_position.is_numeric() &&
           "First column (snp_genomic_position) needs to be decimal with at most 9 digits.");
    assert(c_snp_offset_in_forward_kmer.len <= 2 && c_snp_offset_in_forward_kmer.is_numeric() &&
           "Second column (snp_offset_within_forward_kmer) needs to be decimal with at most 2 digits.");
    assert(c_genome_id.len == 6 && c_genome_id.is_numeric() &&
           "Tenth column (genome_id) needs to be decimal with exactly 6 digits.");
    const auto forward_offset = strtoul(c_snp_offset_in_forward_kmer.c_str, NULL, 10);
    assert(forward_offset >= 0 && forward_offset < K &&
           "Second column (snp_offset_within_forward_kmer) needs to be decimal in range 0 .. K-1, inclusive.");
    for (auto kc : kmer_columns) {
      const char *allele_c_str = kc % 2 ? "1" : "0";     // "0" - major, "1" - minor
      const int reverse_complement = ((kc - 2) / 2) % 2; // 0 - forward, 1 - reverse complement
      const int offset = reverse_complement ? (K - 1 - forward_offset) : forward_offset;
      const string snp_coord_str = string(c_genome_id.c_str) + allele_c_str + c_snp_genomic_position.c_str;
      const uint64_t snp_coord = strtoull(snp_coord_str.c_str(), NULL, 10);
      // 16 decimal digits fit into 56 bits, so we have 8 bits left to encode the rc and snp_offset.
      assert(snp_coord < (1ULL << 56));
      assert(offset < (1ULL << 7));
      // Note that just from the snp_with_rc_and_offset we could recover the kmer if we have the reference genome.
      const uint64_t snp_with_rc_and_offset = (snp_coord << 8) | (reverse_complement << 7) | offset;
      auto kmer = seq_encode<K>(profile_cols[kc].c_str);
      profiles.push_back(KmerData(kmer, snp_with_rc_and_offset));
    }
  }
  fclose(fp);
  sort(profiles.begin(), profiles.end());

  if (profiles.size() == 0) {
    cerr << "Empty file " << f_profile << endl;
    assert(false && "Empty sckmer_profiles file.");
  }

  if (allowed_count == 0) {
    cerr << "Empty file " << f_allowed << endl;
    assert(false && "Empty sckmer_allowed file.");
  }

  // left-join "profiles" into "allowed" results, taking advantage of matching sort order
  auto p = profiles.begin() - 1;
  for (auto a = start; a != end; ++a) {
    const auto akmer = get<0>(*a);
    const auto asnp = get<1>(*a);
    uint64_t pkmer, psnp;
    do {
      ++p;
      assert((p != profiles.end()) && "sckmer_profiles not superset of sckmer_allowed");
      pkmer = get_kmer(*p);
      psnp = get_snp(*p);
      assert(((pkmer < akmer) || ((pkmer == akmer) && (psnp <= asnp))) && "sckmer_profiles not superset of sckmer_allowed");
    } while ((akmer != pkmer) || (asnp != psnp));
    get<1>(*a) = get<1>(*p);
  }

  cerr << (string("Loaded ") + to_string(allowed_count) + " kmers from files " + f_allowed + " and " + f_profile + " in " +
           to_string((chrono_time() - t_start) / 1000.0) + " secs.\n");
}

uint64_t get_species(const KmerData &kd) {
  // We want the most signifficant 6 digits.  The number of digits to begin with is at most 16.
  // Here is a slow way to do it:
  //     x = stoi(to_string(get_snp(kmer_data)).substr(0, 6));
  // The fast way is below.  This actually makes a big difference in overall perf.
  // if it has 13 or more digits, remove the last 7
  auto x = get_snp(kd);
  if (x >= 1000000000000) { // 10**12
    x /= 10000000;          // 10**7
  }
  // if it has 10 or more digits, remove the last 4
  if (x >= 1000000000) { // 10**9
    x /= 10000;          // 10**4
  }
  // so long as it still has 7 or more digits, keep chopping off the last digit
  while (x >= 1000000) { // 10**6
    x /= 10;             // 10**1
  }
  return x;
}

void multi_btc64(int n_path, const char **kpaths) {

  auto timeit = chrono_time();

  assert(n_path > 0);

  // Just enough threads to be able to saturate input I/O, but not too many so
  // the merging of results in memory doesn't kill us.  In practice this produces
  // anywhere between 2x and 4x speedup, depending on the available I/O bandwidth.
  constexpr auto MAX_THREADS = 24;
  int num_running = 0;
  mutex mtx;

  bool running[MAX_THREADS];
  for (auto &r : running) {
    r = false;
  }

  vector<KmerData> kdb;
  vector<KmerData> *kkdb[MAX_THREADS];
  vector<KmerData> *scratch[MAX_THREADS];
  kkdb[0] = &kdb;
  for (int i = 1; i < MAX_THREADS; ++i) {
    kkdb[i] = new vector<KmerData>();
  }
  for (int i = 0; i < MAX_THREADS; ++i) {
    scratch[i] = new vector<KmerData>();
    scratch[i]->resize(3000000);
  }

  // This function runs in each thread.  It loads the bits of an input file,
  // then appends those bits to the kdb vector.
  auto thread_func = [&num_running, &running, &mtx, &kkdb, &scratch](const char *input_file_path, int thread_id) {
    assert(0 <= thread_id && thread_id < MAX_THREADS);
    bit_load(input_file_path, *(kkdb[thread_id]), *(scratch[thread_id]));
    mtx.lock();
    {
      --num_running;
      running[thread_id] = false;
    }
    mtx.unlock();
  };

  // Wait for available thread slot (out of MAX_THREADS), and return its thread_id.
  auto wait_for_thread_slot = [&num_running, &running, &mtx]() -> int {
    bool thread_slot_is_available = false;
    int thread_id = -1;
    while (!(thread_slot_is_available)) {
      mtx.lock();
      {
        thread_slot_is_available = (num_running < MAX_THREADS);
        if (thread_slot_is_available) {
          ++num_running;
          for (int i = 0; i < MAX_THREADS; ++i) {
            if (!(running[i])) {
              running[i] = true;
              thread_id = i;
              break;
            }
          }
        }
      }
      mtx.unlock();
      if (!(thread_slot_is_available)) {
        this_thread::sleep_for(chrono::milliseconds(10));
      }
    }
    assert(0 <= thread_id && thread_id < MAX_THREADS);
    return thread_id;
  };

  for (int i = 0; i < n_path; ++i) {
    auto thread_id = wait_for_thread_slot();
    thread(thread_func, kpaths[i], thread_id).detach();
  }

  // Wait for all threads to complete.
  for (int i = 0; i < MAX_THREADS; ++i) {
    wait_for_thread_slot();
  }

  // Some day we should do this without the merge.  Or, with a merge that's instantaneous.
  // That can be achieved with indirection on the indices of the final sorted array,
  // similar to page table address translation.
  auto time_merge = chrono_time();
  cerr << "Merging input data from " << MAX_THREADS << " threads." << endl;

  for (int i = 0; i < MAX_THREADS; ++i) {
    delete scratch[i];
  }

  // note *kkdb[0] is an alias for kdb
  for (int i = 1; i < MAX_THREADS; ++i) {
    kdb.insert(kdb.end(), kkdb[i]->begin(), kkdb[i]->end());
    delete kkdb[i];
  }

  cerr << "Merged input data.  That took " << (chrono_time() - time_merge) / 1000 << " secs." << endl;
  cerr << "Loaded all files! "
       << "That took " << (chrono_time() - timeit) / 1000 << " secs." << endl;
  cerr << "The unfiltered kmer list has " << kdb.size() << " kmers" << endl;

  cerr << "Sorting..." << endl;
  timeit = chrono_time();
  const auto kdb_size = kdb.size();
  if (kdb_size < 1000 * 1000) {
    sort(kdb.begin(), kdb.end());
  } else {
    // Hacky parallel sort, funny c++ hasn't gotten around to making this standard.
    // Cuts the sort time in half, which is a considerable part of the program's overall runtime.
    using iterator = vector<KmerData>::iterator;
    auto sorter = [](iterator start, iterator end) { sort(start, end); };
    auto merger = [](iterator start, iterator middle, iterator end) {
      // LOL, turns out this "inplace_merge" isn't actually in-place, so it allocs tons of RAM!
      // Should do something much better here --- like parallel quicksort with depth 2.
      inplace_merge(start, middle, end);
    };
    auto q0 = kdb.begin();
    auto q1 = q0 + kdb_size / 4;
    auto q2 = q1 + kdb_size / 4;
    auto q3 = q2 + kdb_size / 4;
    auto q4 = kdb.end();
    thread sorters[] = {
        thread(sorter, q0, q1),
        thread(sorter, q1, q2),
        thread(sorter, q2, q3),
        thread(sorter, q3, q4),
    };
    for (auto &t : sorters) {
      t.join();
    }
    auto m = thread(merger, q0, q1, q2);
    merger(q2, q3, q4);
    m.join();
    merger(q0, q2, q4);
  }
  cerr << "Sorting done! "
       << "It took " << (chrono_time() - timeit) / 1000 << " secs." << endl;

  uint64_t multispecies_kmers = 0;
  uint64_t multispecies_unique_kmers = 0;
  uint64_t monospecies_kmers = 0;
  uint64_t monospecies_unique_kmers = 0;

  cerr << "Starting to check conflicts:  Will filter out kmers that occur in multiple species." << endl;

  timeit = chrono_time();
  ofstream fh(OUT_PATH, ofstream::binary);

  // Scroll through kdb, choosing whether to emit or drop a kmer depending on whether all hits
  // from that kmer are to the same species or to multiple species.
  auto current = kdb.begin();
  while (current != kdb.end()) {
    // Scan forward over all hits from the current kmer value.  Are all those hits to the same species?
    const auto kmer = get_kmer(*current);
    const auto species = get_species(*current);
    bool kmer_is_monospecific = true;
    auto next = current + 1;
    while (next != kdb.end() && get_kmer(*next) == kmer) {
      kmer_is_monospecific = kmer_is_monospecific && (get_species(*next) == species);
      ++next;
    }
    if (kmer_is_monospecific) {
      // Every snp for this kmer is from the same species.  Emit.
      // NOTE:  The snp bytes will come before the kmer bytes in the output.
      // See the comment on the definition of KmerData above at the very top.
      fh.write((char *)&(*current), sizeof(*current) * (next - current));
      monospecies_kmers += (next - current);
      ++monospecies_unique_kmers;
    } else {
      // The current kmer hits multiple species.  Suppress.
      multispecies_kmers += (next - current);
      ++multispecies_unique_kmers;
    }
    current = next;
  }

  cerr << "Filtering done! "
       << "It took " << (chrono_time() - timeit) / 1000 << " secs." << endl;
  cerr << "The filtered kmer list has " << monospecies_kmers << " monospecies kmers (" << monospecies_unique_kmers
       << " unique)." << endl;
  cerr << "Purging conflicts removed " << multispecies_kmers << " multispecies kmers (" << multispecies_unique_kmers
       << " unique)." << endl;

  fh.close();

  if (multispecies_kmers) {
    cerr << "ERROR:  Multispecies kmers found!   There is possibly a problem with the upstream sckmerdb_allowed.tsv "
            "generator.\n";
    assert(false && "multispecies kmers present in input");
  }
}

void display_usage(const char *fname) {
  cout << "usage: " << fname << " fpath [fpath ...]\n"
       << "\n"
       << "where fpath looks like DDDDDD.sckmer_allowed.tsv, for some 6-digit species id DDDDDD, and\n"
       << "a matching file DDDDDD.sckmer_profiles.tsv must exist in the same place for the same DDDDDD.\n";
}

int main(int argc, const char **argv) {
  if ((argc == 2 && (string(argv[1]) == "-h")) || (argc <= 1)) {
    display_usage(argv[0]);
  } else {
    multi_btc64(argc - 1, argv + 1);
  }
  return 0;
}
