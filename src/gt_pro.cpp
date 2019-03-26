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
#include <libgen.h>
#include <regex>

using namespace std;

constexpr auto step_size = 32 * 1024 * 1024;
constexpr auto buffer_size = 32 * 1024 * 1024;

// The DB k-mers are 31-mers.
constexpr auto K = 31;

// 2 bits to encode each ACTG letter
constexpr auto BITS_PER_BASE = 2;

// number of bits to encode entire K-mer
constexpr auto K2 = BITS_PER_BASE * K;

constexpr uint64_t LSB = 1;

// Choose appropriately sized integer types to represent offsets into
// the database.
using LmerRange = uint64_t;
constexpr auto START_BITS = 48;
constexpr auto END_MINUS_START_BITS = 64 - START_BITS;
constexpr auto MAX_START = (LSB << START_BITS) - LSB;
constexpr auto MAX_END_MINUS_START = (LSB << END_MINUS_START_BITS) - LSB;

// This param is only useful for perf testing.  The setting below, not
// to exceed 64 TB of RAM, is equivalent to infinity in 2019.
constexpr auto MAX_MMAP_GB = 64*1024;
constexpr auto MAX_END = MAX_MMAP_GB * (LSB << 30) / 8;

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

template <int M2, int M3>
bool kmer_lookup_work(LmerRange* lmer_indx, uint64_t* mmer_present, uint64_t* data, int channel, char* in_path, char* o_name, int aM2, int aM3){

	if (aM2 != M2 || aM3 != M3) {
		return false;
	}

	const uint64_t MAX_PRESENT = (LSB << M3) - LSB;
    constexpr int XX = (M3 + 1) / 2;  // number DNA letters to cover MAX_PRESENT

	uint8_t code_dict[1 << (sizeof(char) * 8)];
	make_code_dict(code_dict);

	auto out_path = string(o_name) + "." + to_string(channel) + ".tsv";

	//Matching: lmer table lookup then linear search 
	vector<char> buffer(buffer_size);
	char* window = buffer.data();

	uint64_t n_lines = 0;
	
	//  Print progress update every 5 million lines.
	constexpr uint64_t PROGRESS_UPDATE_INTERVAL = 5 * 1000 * 1000;

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

					const auto mmer_pres = seq_encode<uint64_t, XX>(seq_buf + j + K - XX, code_dict) & MAX_PRESENT;
					const auto mpres = mmer_present[mmer_pres / 64] >> (mmer_pres % 64);

					if (mpres & 1) {
						const auto kmer = seq_encode<uint64_t, K>(seq_buf + j, code_dict);
						const auto lmer = kmer >> M2;
						const auto range = lmer_indx[lmer];
						const auto start = range >> END_MINUS_START_BITS;
						const auto end = min(MAX_END, start + (range & MAX_END_MINUS_START));

						for (uint64_t z = start;  z < end;  z += 2) {
							const auto db_kmer = data[z];
							if (kmer == db_kmer) {
								const auto db_snp = data[z + 1];
								if (footprint.find(db_snp) == footprint.end()) {
									kmer_matches.push_back(db_snp);
									footprint.insert({db_snp, 1});
								}
							} else if (kmer < db_kmer) {
								break;
							}
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

	return true;
}

void kmer_lookup(LmerRange* lmer_indx, uint64_t* mmer_present, uint64_t* data, int channel, char* in_path, char* o_name, int M2, const int M3) {
	// Only one of these will really run.  By making them known at compile time, we increase speed.
	// The command line params corresponding to these options are L in {26, 27, 28, 29, 30}  x  M in {30, 32, 34, 35, 36, 37}.
	bool match = false;

	match = match || kmer_lookup_work<32, 30>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<32, 32>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<32, 34>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<32, 35>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<32, 36>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<32, 37>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);

	match = match || kmer_lookup_work<33, 30>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<33, 32>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<33, 34>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<33, 35>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<33, 36>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<33, 37>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);

	match = match || kmer_lookup_work<34, 30>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<34, 32>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<34, 34>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<34, 35>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<34, 36>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<34, 37>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);

	match = match || kmer_lookup_work<35, 30>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<35, 32>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<35, 34>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<35, 35>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<35, 36>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<35, 37>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);

	match = match || kmer_lookup_work<36, 30>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<36, 32>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<36, 34>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<36, 35>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<36, 36>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);
	match = match || kmer_lookup_work<36, 37>(lmer_indx, mmer_present, data, channel, in_path, o_name, M2, M3);

	assert(match && "See comment for supporrted values of L and M.");
}

void display_usage(char* fname){
	cout << "usage: " << fname << " -d <sckmerdb_path: string> -t <n_threads; int; default 1> -o <out_prefix; string; default: cur_dir/out> [-h] input1 [input2 ...]\n";
}


template <class ElementType>
struct DBIndex {

	std::vector<ElementType> elements;
	uint64_t expected_element_count;
	std::string filename;
	bool loaded;

	DBIndex(const string& filename, const uint64_t expected_element_count=0)
		: filename(filename),
		  loaded(false),
		  expected_element_count(expected_element_count)
	{}

	ElementType* address() {
		assert(elements.size() > 0);
		return &(elements[0]);
	}

	void allocate_then_load() {
		assert(!(loaded));
		FILE* dbin = fopen(filename.c_str(), "rb");
		if (dbin) {
			cerr << chrono_time() << ":  Loading " << filename << endl;
			size_t filesize;
			if (expected_element_count) {
				filesize = expected_element_count * sizeof(ElementType);
			} else {
				filesize = get_fsize(filename.c_str());
				expected_element_count = filesize / sizeof(ElementType);
			}
			elements.resize(expected_element_count);
			const auto loaded_element_count = fread(address(), sizeof(ElementType), expected_element_count, dbin);
			if (loaded_element_count == expected_element_count) {
				cerr << chrono_time() << ":  Loaded " << filename << endl;
				loaded = true;
			} else {
				cerr << chrono_time() << ":  Failed to load " << filename << ".  This is fine, it will only slow down init." << endl;
			}
			fclose(dbin);
		}
		if (!(loaded)) {
			cerr << chrono_time() << ": " << filename << " will be recomputed." << endl;
			elements.resize(expected_element_count);
		}
	}

	void save_if_did_not_exist() {
		if (!(loaded)) {
			auto l_start = chrono_time();
			FILE* dbout = fopen(filename.c_str(), "wb");
			assert(dbout);
			const auto saved_element_count = fwrite(address(), sizeof(ElementType), elements.size(), dbout);
			fclose(dbout);
			assert(saved_element_count == elements.size());
			cerr << chrono_time() << ":  Done writing " << filename << ". That took " << (chrono_time() - l_start) / 1000 << " more seconds." << endl;
		}
	}
};


int main(int argc, char** argv) {
	extern char *optarg; 
	extern int optind;

	bool dbflag = false;
	bool inflag = false;

	char* fname = argv[0];
	char* db_path = (char *)"";
	char* oname = (char *)"./out";

	// Number of bits in the prefix part of the K-mer (also called L-mer,
	// even though it might not correspond to an exact number of bases).
	// This has a substantial effect on memory use.  Rule of thumb for
	// perf is L2 >= K2 - M3.	
	auto L2 = 29;

	// Number of bits in the MMER_PRESENT index.  This has a substantial effect
	// on memory use.  Rule of thumb for perf is M3 >= 4 + log2(DB cardinality).
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
			case 'h': case '?':
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
	//assert(L2 >= K2 - M3);

	const auto LMER_MASK = (LSB << L2) - LSB;
	const auto MMER_MASK = (LSB << M2) - LSB;
	const auto MAX_PRESENT = (LSB << M3) - LSB;

	cout << fname << '\t' << db_path << '\t' << n_threads << "\t" << (preload ? "preload" : "mmap") << "\t" << L2 << "\t" << M3 << endl;

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
	cerr << chrono_time() << ":  " << "Starting to load DB: " << db_path << endl;

	size_t filesize = get_fsize(db_path);
	//Open file
	int fd = open(db_path, O_RDONLY, 0);
	assert(fd != -1);
	uint64_t* mmappedData = (uint64_t *) mmap(NULL, filesize, PROT_READ, MMAP_FLAGS, fd, 0);
	assert(mmappedData != MAP_FAILED);

	string dbbase = string(basename(db_path));
	dbbase = regex_replace(dbbase, regex("\\.bin$"), "");
	dbbase = regex_replace(dbbase, regex("\\."), "_");

	const auto OPTIMIZED_DB_MMER_PRESENT = dbbase + "_optimized_db_mmer_present_" + to_string(M3) + ".bin";
    DBIndex<uint64_t> db_mmer_present(OPTIMIZED_DB_MMER_PRESENT, (1 + MAX_PRESENT) / 64);  // 1 bit per possible mmer
	db_mmer_present.allocate_then_load();
	const bool recompute_mmer_present = !(db_mmer_present.loaded);
	uint64_t* mmer_present = db_mmer_present.address();

	const auto OPTIMIZED_DB_LMER_INDX = dbbase + "_optimized_db_lmer_indx_" + to_string(L2) + ".bin";
    DBIndex<LmerRange> db_lmer_index(OPTIMIZED_DB_LMER_INDX, 1 + LMER_MASK);
	db_lmer_index.allocate_then_load();
	const bool recompute_lmer_indx = !(db_lmer_index.loaded);
	LmerRange* lmer_indx = db_lmer_index.address();

	uint64_t last_lmer;
	uint64_t start = 0;

	uint64_t lmer_count = -1;

	auto data = mmappedData;
	if (preload) {
		cerr << chrono_time() << ":  Preloading entire DB as requested.  Usually this makes performance worse." << endl;
		data = new uint64_t[filesize / 8];
		for (uint64_t i = 0;  i < filesize / 8;  ++i) {
			data[i] = mmappedData[i];
		}
		cerr << chrono_time() << ":  Preloaded entire DB as requested." << endl;
	}

	if (recompute_mmer_present || recompute_lmer_indx) {
		lmer_count = filesize ? 1 : 0;
		for (uint64_t end = 0;  end < filesize / 8;  end += 2) {
			const auto kmer = data[end];
			const auto lmer = kmer >> M2;
			if (recompute_mmer_present) {
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
			if (recompute_lmer_indx) {
				lmer_indx[lmer] = (start << END_MINUS_START_BITS) | (end - start + 1);
			}
			last_lmer = lmer;
		}
	}

	db_lmer_index.save_if_did_not_exist();
	db_mmer_present.save_if_did_not_exist();

	cerr << chrono_time() << ":  " << "Done with init for DB with " << (filesize / 16) << " mmers.  That took " << (chrono_time() - l_start) / 1000 << " seconds." << endl;

	l_start = chrono_time();

	vector<thread> th_array;
	int tmp_counter = 0;
	for(; optind < argc; optind++) {
		th_array.push_back(thread(kmer_lookup, lmer_indx, mmer_present, data, optind - in_pos, argv[optind], oname, M2, M3));
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

	if (preload) {
		delete [] data;
	}

	cerr << chrono_time() << ":  " << " Totally done: " << (chrono_time() - l_start) / 1000 << " seconds elapsed processing reads, after DB was loaded."  << endl;

	return 0;
}
