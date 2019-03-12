#if __linux__
#include <linux/version.h>
#if LINUX_VERSION_CODE > KERNEL_VERSION(2,6,22)
#define _MAP_POPULATE_AVAILABLE
#endif
#endif

#ifdef _MAP_POPULATE_AVAILABLE
#define MMAP_FLAGS (MAP_PRIVATE | MAP_POPULATE)
#else
#define MMAP_FLAGS MAP_PRIVATE
#endif

#include <sys/mman.h>
#include <stdio.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

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

using namespace std;

constexpr auto step_size = 256 * 1024 * 1024;
constexpr auto buffer_size = 256 * 1024 * 1024;

// bits per base
constexpr int bpb = 2;

size_t get_fsize(const char* filename) {
	struct stat st;
	stat(filename, &st);
	return st.st_size;
}

template <class int_type>
int_type bit_encode(const char c) {
	switch (c) {
		case 'A': case 'a': return 0;
		case 'C': case 'c': return 1;
		case 'G': case 'g': return 2;
		case 'T': case 't': return 3;
	}

	assert(false);
}

template <class int_type>
char bit_decode(const int_type bit_code) {
	switch (bit_code) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
	}
	assert(false);
}

template <class int_type>
void make_code_dict(int_type* code_dict) {
	code_dict['A'] = bit_encode<int_type>('A');
	code_dict['a'] = bit_encode<int_type>('a');
	code_dict['C'] = bit_encode<int_type>('C');
	code_dict['c'] = bit_encode<int_type>('c');
	code_dict['G'] = bit_encode<int_type>('G');
	code_dict['g'] = bit_encode<int_type>('g');
	code_dict['T'] = bit_encode<int_type>('T');
	code_dict['t'] = bit_encode<int_type>('t');
}

template <class int_type>
int_type seq_encode(const char* buf, int len, const int_type* code_dict, const int_type b_mask) {
	int_type seq_code = 0;
	for (int i=0;  i < len;  ++i) {
		const int_type b_code = code_dict[buf[i]];
		seq_code |= ((b_code & b_mask) << (bpb * (len - i - 1)));
	}
	return seq_code;
}

template <class int_type>
void seq_decode(char* buf, const int len, const int_type seq_code, const int_type b_mask) {
	for (int i=0;  i < len-1;  ++i) {
		const int_type b_code = (seq_code >> (bpb * (len - i - 2))) & b_mask;
		buf[i] = bit_decode<int_type>(b_code);
	}

	buf[len-1] = '\0';
}

long chrono_time() {
	using namespace chrono;
	return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

using LmerRange = tuple<uint64_t, uint64_t>;

void kmer_lookup(LmerRange* lmer_indx, vector<uint32_t>& mmers, vector<uint64_t>& snps, int channel, char* in_path, char* o_name){
	uint32_t lsb = 1;
	uint32_t b_mask = (lsb << bpb) - lsb;

	uint32_t code_dict[1 << (sizeof(char) * 8)];
	make_code_dict<uint32_t>(code_dict);

	auto out_path = string(o_name) + "." + to_string(channel) + ".tsv";

	constexpr int L = 15;
	constexpr int K = 31;
	constexpr int M = K - L;

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
	char lmer_buf[L + 1];
	char mmer_buf[M + 1];

	lmer_buf[L] = '\0';
	mmer_buf[L] = '\0';

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

					// lmer_buf[0:L] := seq_buf[j:j+L]
					for (int z = j;  z < j + L;  ++z) {
						lmer_buf[z - j] = seq_buf[z];
					}

					// ensured by prior initialization
					assert(lmer_buf[L] == '\0');

					const auto lcode = seq_encode<uint32_t>(lmer_buf, L, code_dict, b_mask);
					const auto range = lmer_indx[lcode];
					const auto start = get<0>(range);
					const auto end = get<1>(range);

					if (end) {

						// mmer_buf[0:K-L] := seq_buf[j+L:j+K]
						for (int z = j + L;  z < j + K;  ++z) {
							mmer_buf[z - j -L] = seq_buf[z];
						}

						// ensured by prior initialization
						assert(mmer_buf[M] == '\0');

						auto mcode = seq_encode<uint32_t>(mmer_buf, M, code_dict, b_mask);

						// linear search
						for (uint64_t z = start; z < end; ++z) {
							if (mcode == mmers[z]) {
								if (footprint.find(snps[z]) != footprint.end()) {
									//do nothing
								} else {
									kmer_matches.push_back(snps[z]);
									footprint.insert({snps[z], 1});
								}
							} else {
								// cout << "    no match: " << mcode << " - " << mmers[j] << '\n';
							}
						}
					}
				}
			}

			// next token, please
			footprint.clear();
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

	uint64_t start = 0;
	uint64_t end = 0;

	vector<uint32_t> mmers;
	vector<uint64_t> snps;
	uint32_t last_lmer = 2147483648;  // 1u << 31u, this value doesn't matter

	mmers.reserve(filesize / 8);
	snps.reserve(filesize / 8);

	auto l_start = chrono_time();
	cerr << chrono_time() << ":  " << "Starting to load DB: " << db_path << endl;

	constexpr uint64_t BILLION = ((uint64_t) 1) << (uint64_t) 30;  // 2 ** 30
	LmerRange *lmer_indx = new LmerRange[BILLION]();

	cerr << chrono_time() << ":  " << "Allocated memory.  That took " << (chrono_time() - l_start) / 1000 << " seconds." << endl;
	l_start = chrono_time();

	for (uint64_t i = 0;  i < filesize / 8;  i += 2) {
		auto kmer = mmappedData[i];
		uint32_t lmer = (uint32_t)((kmer & 0x3FFFFFFF00000000LL) >> 32);
		uint32_t mmer = (uint32_t)(kmer & 0xFFFFFFFFLL);
		mmers.push_back(mmer);
		++end; 
		snps.push_back(mmappedData[i+1]);
		if (i > 0 && lmer != last_lmer) {
			start = end - 1;
		}
		// Invariant:  The data loaded so far for lmer reside at offsets start, start+1, ..., end-1.
		lmer_indx[lmer] = make_tuple(start, end);
		last_lmer = lmer;
	}

	cerr << chrono_time() << ":  " << "Done loading DB.  That took " << (chrono_time() - l_start) / 1000 << " more seconds." << endl;
	l_start = chrono_time();

	int rc = munmap(mmappedData, filesize);
	assert(rc == 0);
	close(fd);
	

	vector<thread> th_array;
	int tmp_counter = 0;
	for(; optind < argc; optind++) {
		th_array.push_back(thread(kmer_lookup, lmer_indx, ref(mmers), ref(snps), optind - in_pos, argv[optind], oname));
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

	cerr << chrono_time() << ":  " << " Totally done: " << (chrono_time() - l_start) / 1000 << " seconds elapsed processing reads, after DB was loaded."  << endl;

	return 0;
}
