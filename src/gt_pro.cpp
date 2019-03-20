#if __linux__
#include <linux/version.h>
#if LINUX_VERSION_CODE > KERNEL_VERSION(2,6,22)
#define _MAP_POPULATE_AVAILABLE
#endif
#endif

#ifdef _MAP_POPULATE_AVAILABLE
#define MMAP_FLAGS (MAP_PRIVATE | MAP_POPULATE)
#else
#define MMAP_FLAGS (MAP_PRIVATE)
#endif

#include <sys/mman.h>
#include <stdio.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>
#include <string.h>

#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <thread>
#include <limits>

using namespace std;

constexpr auto step_size = 256 * 1024 * 1024;
constexpr auto buffer_size = 256 * 1024 * 1024;

// See bit_encode
constexpr uint8_t BITS_PER_BASE = 2;

constexpr int L = 15;
constexpr int K = 31;
constexpr int M = K - L;

static_assert(L > 0);
static_assert(M > 0);

constexpr uint64_t LMER_MASK = (((uint64_t) 1) << (BITS_PER_BASE * L)) - 1;
constexpr uint64_t MMER_MASK = (((uint64_t) 1) << (BITS_PER_BASE * M)) - 1;

// Given L and M above, choose appropriately-sized integer types
// to represent lmers and mmers.
using LmerType = uint32_t;
using MmerType = uint32_t;

static_assert(LMER_MASK <= numeric_limits<LmerType>::max());
static_assert(MMER_MASK <= numeric_limits<MmerType>::max());

const uint64_t MAX_PRESENT = (((uint64_t) MMER_MASK) << 4) | 0xf;

// Choose appropriately sized integer types to represent offsets into
// the array of all mmers.
using LmerRange = uint64_t;
constexpr auto MAX_START = (((uint64_t) 1) << 48) - 1;
constexpr auto MAX_END_MINUS_START = 65535;

size_t get_fsize(const char* filename) {
	struct stat st;
	stat(filename, &st);
	return st.st_size;
}

uint8_t bit_encode(const char c) {
	switch (c) {
		case 'A': case 'a': return 0;
		case 'C': case 'c': return 1;
		case 'G': case 'g': return 2;
		case 'T': case 't': return 3;
	}
	assert(false);
}

char bit_decode(const uint8_t bit_code) {
	switch (bit_code) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
	}
	assert(false);
}

void make_code_dict(uint8_t* code_dict) {
	code_dict['A'] = bit_encode('A');
	code_dict['a'] = bit_encode('a');
	code_dict['C'] = bit_encode('C');
	code_dict['c'] = bit_encode('c');
	code_dict['G'] = bit_encode('G');
	code_dict['g'] = bit_encode('g');
	code_dict['T'] = bit_encode('T');
	code_dict['t'] = bit_encode('t');
}

template <class int_type, int len>
int_type seq_encode(const char* buf, const uint8_t* code_dict) {
	int_type seq_code = 0;
	// This loop may be unrolled by the compiler because len is a compile-time constant.
	for (int bitpos = (len - 1) * BITS_PER_BASE;  bitpos >= 0;  bitpos -= BITS_PER_BASE) {
		const int_type b_code = code_dict[*buf++];
		seq_code |= (b_code << bitpos);
	}
	return seq_code;
}

template <class int_type>
void seq_decode(char* buf, const int len, const int_type seq_code) {
	constexpr uint8_t B_MASK = (1 << BITS_PER_BASE) - 1;
	for (int i=0;  i < len - 1;  ++i) {
		const uint8_t b_code = B_MASK & (seq_code >> (BITS_PER_BASE * (len - i - 2)));
		buf[i] = bit_decode(b_code);
	}
	buf[len - 1] = '\0';
}

long chrono_time() {
	using namespace chrono;
	return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

void kmer_lookup(LmerRange* lmer_indx, uint64_t* mmer_present, uint64_t* mmappedData, int channel, char* in_path, char* o_name){
	uint8_t code_dict[1 << (sizeof(char) * 8)];
	make_code_dict(code_dict);

	auto out_path = string(o_name) + "." + to_string(channel) + ".tsv";

	//Matching: lmer table lookup then linear search 
	vector<char> buffer(buffer_size);
	char* window = buffer.data();

	uint64_t n_lines = 0;
	
	// Print progress update every 5 million lines.
	constexpr uint64_t PROGRESS_UPDATE_INTERVAL = 5*1000*1000;

	// Reads that contain wildcard characters ('N' or 'n') are split into
	// tokens at those wildcard characters.  Each token is processed as
	// though it were a separate read.
	constexpr int MAX_TOKEN_LENGTH = 500;
	constexpr int MIN_TOKEN_LENGTH = 31;

	char seq_buf[MAX_TOKEN_LENGTH];

	// This ranges from 0 to the length of the longest read (could exceed MAX_TOKEN_LENGTH).
	int token_length = 0;

	vector<uint64_t> kmer_matches;

	unordered_map<uint64_t, int> footprint;

	int fd = open(in_path, O_RDONLY);

	auto s_start = chrono_time();
	char c = '\0';

	while (true) {
		const ssize_t bytes_read = read(fd, window, step_size);

		if (bytes_read == 0)
			break;

		if (bytes_read == (ssize_t) - 1) {
			cerr << chrono_time() << ":  " << "unknown fatal error, when read stdin input" << endl;
			exit(EXIT_FAILURE);
		}

		for (uint64_t i = 0;  i < bytes_read;  ++i) {

			if (c == '\n') {
				++n_lines;
				if ((n_lines + 1) % PROGRESS_UPDATE_INTERVAL == 0) {
					cerr << chrono_time() << ":  " << ((n_lines + 3) / 4) << " reads were scanned after "
						 << (chrono_time() - s_start) / 1000 << " seconds from file "
						 << in_path << endl;
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
				for (int j = 0;  j <= token_length - K;  ++j) {

					const uint64_t kmer = seq_encode<uint64_t, K>(seq_buf + j, code_dict);
					const uint64_t mmer_pres = kmer & MAX_PRESENT;
					const auto mpres = mmer_present[mmer_pres / 64] >> (mmer_pres % 64);

					if (mpres & 1) {
						const LmerType lmer = kmer >> (M * BITS_PER_BASE);
						const auto range = lmer_indx[lmer];
						const auto start = range >> 16;
						const auto end = start + (range & 0xFFFF);

						for (uint64_t z = start;  z < end;  z += 2) {
							const auto db_kmer = mmappedData[z];
							if (kmer == db_kmer) {
								const auto db_snp = mmappedData[z + 1];
								if (footprint.find(db_snp) == footprint.end()) {
									kmer_matches.push_back(db_snp);
									footprint.insert({db_snp, 1});
								}
							} //else if (kmer < db_kmer) {
							//	break;
							//}
						}
					}
				}
			}

			// clear footprint for every read instead of every token
			if(c == '\n') {
				footprint.clear();
			}

			// next token, please
			token_length = 0;
		}
	}

	if (token_length != 0) {
		cerr << chrono_time() << ":  " << "Error:  Truncated read sequence at end of file: " << in_path << endl;
		exit(EXIT_FAILURE);
	}

	cerr << chrono_time() << ":  " << "[Done] searching is completed, emitting results for " << in_path << endl;
	ofstream fh(out_path, ofstream::out | ofstream::binary);

	if (kmer_matches.size() == 0) {
		cerr << chrono_time() << ":  " << "zero hits" << endl;	
	} else {
		sort(kmer_matches.begin(), kmer_matches.end());

		uint64_t cur_count = 0;
		uint64_t cur_snp = kmer_matches[0];
		for(auto it = kmer_matches.begin(); it != kmer_matches.end(); it=it+1){
			if (*it != cur_snp) {
				fh << cur_snp << '\t' << cur_count << '\n';
				cur_snp = *it;
				cur_count = 1;
			} else {
				++cur_count;
			}

			// seq_decode(dcd_buf, k, it->first, code_dict, b_mask);
			// fh << *it << "\t" << *(it+1) << "\n";
		}
	}
	cerr << chrono_time() << ":  " << "Completed output for " << in_path << endl;

	fh.close();
}

void display_usage(char* fname){
	cout << "usage: " << fname << " -d <sckmerdb_path: string> -t <n_threads; int; default 1> -o <out_prefix; string; default: cur_dir/out> [-h] input1 [input2 ...]\n";
}

int main(int argc, char** argv) {
	extern char *optarg; 
	extern int optind;

	bool dbflag = false;
	bool inflag = false;

	char* fname = argv[0];
	char* db_path = (char *)"";
	char* oname = (char *)"./out";

	int n_threads = 1;

	int opt;
	while ((opt = getopt(argc, argv, "d:t:o:h")) != -1) {
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
			case 'h': case '?':
				display_usage(fname);
				exit(1);
		}	
	}

	cout << fname << '\t' << db_path << '\t' << n_threads << endl;

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


	size_t filesize = get_fsize(db_path);
	//Open file
	int fd = open(db_path, O_RDONLY, 0);
	assert(fd != -1);
	uint64_t* mmappedData = (uint64_t *) mmap(NULL, filesize, PROT_READ, MMAP_FLAGS, fd, 0);
	assert(mmappedData != MAP_FAILED);

	auto l_start = chrono_time();
	cerr << chrono_time() << ":  " << "Starting to load DB: " << db_path << endl;

	LmerRange *lmer_indx = new LmerRange[1 + LMER_MASK];
	uint64_t *mmer_present = new uint64_t[(1 + MAX_PRESENT) / 64];  // allocate 1 bit per possible mmer

	const char* OPTIMIZED_DB_MMER_PRESENT = "optimized_db_mmer_present.bin";
	const char* OPTIMIZED_DB_LMER_INDX = "optimized_db_lmer_indx.bin";

	FILE* dbin = fopen(OPTIMIZED_DB_MMER_PRESENT, "rb");
	if (dbin) {
		int success = fread(mmer_present, 8, (1 + MAX_PRESENT) / 64, dbin);
		if (success == (1 + MAX_PRESENT) / 64) {
			cerr << chrono_time() << ":  Read MMER_PRESENT from " << OPTIMIZED_DB_MMER_PRESENT << endl;
			fclose(dbin);
		} else {
			cerr << chrono_time() << ":  Failed to read " << OPTIMIZED_DB_MMER_PRESENT << ".  This is fine, it will only slow down init." << endl;
			fclose(dbin);
			dbin = NULL;
		}
	}

	if (!(dbin)) {
		cerr << chrono_time() << ":  MMER_PRESENT will be recomputed." << endl;
		memset(mmer_present, 0, (1 + MAX_PRESENT) / 8);
	}

	FILE* lmerdbin = fopen(OPTIMIZED_DB_LMER_INDX, "rb");
	if (lmerdbin) {
		int success = fread(lmer_indx, sizeof(LmerRange), (1 + LMER_MASK), lmerdbin);
		if (success == (1 + LMER_MASK)) {
			cerr << chrono_time() << ":  Read LMER_INDX from " << OPTIMIZED_DB_LMER_INDX << endl;
			fclose(lmerdbin);
		} else {
			cerr << chrono_time() << ":  Failed to read " << OPTIMIZED_DB_LMER_INDX << ".  This is fine, it will only slow down init." << endl;
			fclose(lmerdbin);
			lmerdbin = NULL;
		}
	}

	if (!(lmerdbin)) {
		cerr << chrono_time() << ":  LMER_INDX will be recomputed." << endl;
		memset(lmer_indx, 0, sizeof(LmerRange) * (1 + LMER_MASK));
	}

	LmerType last_lmer;
	uint64_t start = 0;

	uint64_t lmer_count = -1;

	if (!(lmerdbin) || !(dbin)) {
		lmer_count = filesize ? 1 : 0;
		for (uint64_t end = 0;  end < filesize / 8;  end += 2) {
			auto kmer = mmappedData[end];
			MmerType mmer = kmer & MMER_MASK;
			LmerType lmer = kmer >> (M * BITS_PER_BASE);
			if (!(dbin)) {
				uint64_t mpres = kmer & MAX_PRESENT;
				mmer_present[mpres / 64] |= ((uint64_t) 1) << (mpres % 64);
			}
			if (end > 0 && lmer != last_lmer) {
				start = end;
				++lmer_count;
			}
			// Invariant:  The data loaded so far for lmer reside at offsets start, start+1, ..., end-1.
			assert(start <= MAX_START);
			assert(end - start < MAX_END_MINUS_START);
			assert(lmer <= LMER_MASK);
			if (!(lmerdbin)) {
				*((uint64_t*)(&(lmer_indx[lmer]))) = (start << 16) | (end - start + 1);
			}
			last_lmer = lmer;
		}
	}

	cerr << chrono_time() << ":  " << "Done loading DB of " << (filesize / 16) << " mmers.  That took " << (chrono_time() - l_start) / 1000 << " seconds." << endl;

	if (!(dbin)) {
		l_start = chrono_time();
		FILE* dbout = fopen(OPTIMIZED_DB_MMER_PRESENT, "wb");
		assert(dbout);
		int success = fwrite(mmer_present, 8, (1 + MAX_PRESENT) / 64, dbout);	
		fclose(dbout);
		assert(success == (1 + MAX_PRESENT) / 64);
		cerr << chrono_time() << ":  Done writing optimized MMER_PRESENT.  That took " << (chrono_time() - l_start) / 1000 << " more seconds." << endl;
	}

	if (!(lmerdbin)) {
		l_start = chrono_time();		
		FILE* dbout = fopen(OPTIMIZED_DB_LMER_INDX, "wb");
		assert(dbout);
		int success = fwrite(lmer_indx, sizeof(LmerRange), (1 + LMER_MASK), dbout);	
		fclose(dbout);
		assert(success == (1 + LMER_MASK));
		cerr << chrono_time() << ":  Done writing optimized LMER_INDX with " << lmer_count << " lmers.  That took " << (chrono_time() - l_start) / 1000 << " more seconds." << endl;
	}

	l_start = chrono_time();

	vector<thread> th_array;
	int tmp_counter = 0;
	for(; optind < argc; optind++) {
		th_array.push_back(thread(kmer_lookup, lmer_indx, mmer_present, mmappedData, optind - in_pos, argv[optind], oname));
		++tmp_counter;

		if (tmp_counter >= n_threads) {

			cerr << chrono_time() << ":  " << "Waiting on all threads from this round to finish before dispatching next round." << endl;

			for (thread & ith : th_array) {
				ith.join();
				// cerr << chrono_time() << ":  " << " Thread joined " << endl;
			}
			
			cerr << chrono_time() << ":  " << "Ready to dispatch next round of threads." << endl;

			th_array.clear();
			tmp_counter = 0;
		} 
	}

	if (th_array.size() > 0){
		for (thread & ith : th_array) {
			ith.join();
			// cerr << chrono_time() << ":  " << " Thread joined " << endl;
		}
		th_array.clear();
	}

	int rc = munmap(mmappedData, filesize);
	assert(rc == 0);
	close(fd);

	cerr << chrono_time() << ":  " << " Totally done: " << (chrono_time() - l_start) / 1000 << " seconds elapsed processing reads, after DB was loaded."  << endl;

	return 0;
}
