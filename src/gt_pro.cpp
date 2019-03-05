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

void kmer_lookup(tuple<uint64_t, uint64_t>* lmer_indx, vector<uint32_t>& mmers, vector<uint64_t>& snps, int channel, char* in_path, char* o_name){
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

	uintmax_t n_lines = 0;
	uintmax_t n_pause = 0;

	int cur_pos = 0;

	constexpr int MAX_READ_LENGTH = 500;
	constexpr int MIN_READ_LENGTH = 50;

	char seq_buf[MAX_READ_LENGTH];
	char lmer_buf[L + 1];
	char mmer_buf[M + 1];

	lmer_buf[L] = '\0';
	mmer_buf[L] = '\0';

	int l_label = 2;
	bool has_wildcard = false;
	bool irregular_read = false;

	vector<uint64_t> kmer_matches;

	unordered_map<uint64_t, int> footprint;

	int fd = open(in_path, O_RDONLY);

	auto s_start = chrono_time();
	while (true) {
		const ssize_t bytes_read = read(fd, window, step_size);

		if (bytes_read == 0)
			break;

		if (bytes_read == (ssize_t) - 1) {
			cerr << chrono_time() << ":  " << "unknown fatal error, when read stdin input" << endl;
			exit(EXIT_FAILURE);
		}

		for (uint64_t i = 0;  i < bytes_read;  ++i) {
			const char c = window[i];

			// In the FASTQ format, every 4 lines define a read.  The first line is the
			// read header.  The next line is the read sequence.  We only care about
			// the read sequences.  We track this by incrementing l_label modulo 4 on
            // each line, then ignore lines for which l_label != 3.  We could just
            // as easily use the equivalent condition (n_lines % 4 != 1).
			assert(l_label == (n_lines + 2) % 4);

            if (c != '\n') {
				// does the current line contain a read sequence?
				if (l_label == 3) {
					// yes, buffer the read's characters into seq_buf[0...MAX_READ_LENGTH-1]
					if (c == 'N') {
						has_wildcard = true;
					}
					if (cur_pos < MAX_READ_LENGTH) {
						seq_buf[cur_pos++] = c;
					} else {
						irregular_read == true;
					}
				}
				// next character, please
				continue;
			}

			assert(c == '\n');

			++n_lines;
			++n_pause;

			if (n_pause == 5*1000*1000) {
				cerr << chrono_time() << ":  " << n_lines << " lines were scanned after "<< (chrono_time() - s_start) / 1000 << " seconds from file " << in_path << endl;
				n_pause = 0;
			}

			if (l_label == 3) {

				assert(cur_pos <= MAX_READ_LENGTH);

				if (cur_pos < MIN_READ_LENGTH) {
					irregular_read = true;	
				}

				if (has_wildcard || irregular_read) {
					has_wildcard = false;
					irregular_read = false;
					l_label = cur_pos = 0;
					footprint.clear();
					continue;
				}

				for (int j = 0;  j <= cur_pos - K;  ++j) {

					// lmer_buf[0:L] := seq_buf[j:j+L]
					for (int z = j;  z < j + L;  ++z) {
						lmer_buf[z - j] = seq_buf[z];
					}

					// ensured by prior initialization
					assert(lmer_buf[L] == '\0');

					auto lcode = seq_encode<uint32_t>(lmer_buf, L, code_dict, b_mask);
					if (get<1>(lmer_indx[lcode])) {

						// mmer_buf[0:K-L] := seq_buf[j+L:j+K]
						for (int z = j + L;  z < j + K;  ++z) {
							mmer_buf[z - j -L] = seq_buf[z];
						}

						// ensured by prior initialization
						assert(mmer_buf[M] == '\0');

						auto mcode = seq_encode<uint32_t>(mmer_buf, M, code_dict, b_mask);

						auto coord = lmer_indx[lcode];
						auto start = get<0>(coord);
						auto end = get<1>(coord);


						/* surprisingly slower
						// binary search
						while (start < end) {
							uint64_t mid = start + (end - start) / 2;

							// cerr << chrono_time() << ":  " << start << '\t' << end << '\t' << mid << '\n';
							if (mmers[mid] == mcode) {
								if (footprint.find(snps[mid]) != footprint.end()) {
									//do nothing
								} else {
									kmer_matches.push_back(snps[mid]);
									footprint.insert({snps[mid], 1});
								}

								break;
							}

							if (mmers[mid] < mcode){
								start = mid + 1;
							} else {
								end = mid;
							}
						}
						*/

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

				footprint.clear();
				l_label = cur_pos = 0;
			} else {
				++l_label;
			}
		}
	}

	if (cur_pos != 0) {
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
	//Execute mmap
	//uint64_t* mmappedData = (uint64_t *) mmap(NULL, filesize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);
	uint64_t* mmappedData = (uint64_t *) mmap(NULL, filesize, PROT_READ, MMAP_FLAGS, fd, 0);
	assert(mmappedData != MAP_FAILED);
	//Write the mmapped data to stdout (= FD #1)

	// write(1, mmappedData, filesize);

	// char seq_buf[k+1];

	uint64_t start = -1;
	uint64_t end = 0;  // FIXME: this init breaks the loop invariant;  should be -1

	vector<uint32_t> mmers;
	vector<uint64_t> snps;
	uint32_t cur_lmer = 2147483648;  // 1u << 31u

	auto l_start = chrono_time();

	cerr << chrono_time() << ":  " << "[OK] Reserving memory." << endl;
	mmers.reserve(filesize / 8);
	snps.reserve(filesize / 8);

	cerr << chrono_time() << ":  " << "[OK] Zeroing index." << endl;

	// unordered_map<uint32_t, tuple<uint64_t, uint64_t>> lmer_indx;
	constexpr uint64_t k_1_BILLION = ((uint64_t) 1) << (uint64_t) 30;
	tuple<uint64_t, uint64_t> *lmer_indx = new tuple<uint64_t, uint64_t>[k_1_BILLION]();

	cerr << chrono_time() << ":  " << "[OK] start to load DB: " << db_path << endl;

	for (uint64_t i = 0;  i < filesize / 8;  i += 2) {
		// seq_decode<uint_fast64_t>(seq_buf, k, mmappedData[i], b_mask);

		auto kmer_int = mmappedData[i];

		uint32_t lmer_int = (uint32_t)((kmer_int & 0x3FFFFFFF00000000LL) >> 32);
		uint32_t mmer_int = (uint32_t)(kmer_int & 0xFFFFFFFFLL);

		mmers.push_back(mmer_int);
		++end;

		snps.push_back(mmappedData[i+1]);

		// cerr << chrono_time() << ":  " << lmer << " - " << mmappedData[i] << '\n';

		if (cur_lmer != lmer_int) {
			if (start == -1){
				start = 0;
			} else {
				// cerr << chrono_time() << ":  " << cur_lmer << ": (" << start << " , " << end << "}\n";
				// cerr << chrono_time() << ":  " << end - start << '\n';
				// mmer/snps entries at position start, start+1, ..., end-1 belong to cur_lmer
				// mmer/snps entry at position end belongs to the next lmer_int
			    auto coord = make_tuple(start, end);
				lmer_indx[cur_lmer] = coord;

				// start, end := mmers.size() - 1, snps.size() - 1
				// now they point to the last element of mmers and snps
				// so the complicated assignents below should really just be "start = end"
				// assert(end == i / 2);
				start = i/2;
				end = start;
			}

			cur_lmer = lmer_int;
		}
	}

    // ^^^ Did we just completely drop the last lmer_int?

	cerr << chrono_time() << ":  " << "Done loading DB.  That took " << (chrono_time() - l_start) / 1000 << " seconds." << endl;
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
