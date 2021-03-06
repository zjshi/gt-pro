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
#include <sys/stat.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <cstring>

#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

using namespace std;


// this simple program scans its binary input (64-bit integers) for k-mers and snp IDs, 
// then query each k-mer in the k-mer pools. 
// Any k-mer found in a k-mer pool is excluded.
// input, output and k-mer pool share the same format: 
// interger1[64-bit; k-mer 1],interger2[64-bit; snp ID 1],interger3[64-bit; k-mer 2],interger4[64-bit; snp ID 2]...

// usage:
//    g++ -O3 --std=c++11 -o db_uniq db_uniq.cpp
//    ./db_uniq -d <binary input> -o <binary output> -L <a list of paths of k-mer pools>
//

// global variable declaration starts here
constexpr auto k = 31;

// set operation mode
// valid values: 0, 1, 2
// 0 is set union operation; 1 is set intersection operation; 2 is set difference([set1-set2]);
constexpr auto s_mod = 0;

// parameters for <unistd.h> file read; from the source of GNU coreutils wc
constexpr auto step_size = 256 * 1024 * 1024;
constexpr auto buffer_size = 256 * 1024 * 1024;


// get time elapsed since when it all began in milliseconds.
long chrono_time() {
    using namespace chrono;
    return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

// number of bits per single nucleotide base
constexpr int bpb = 2;

size_t get_fsize(const char* filename) {
	struct stat st;
	stat(filename, &st);
	return st.st_size;
}

char* get_ftype(const char* filename) {
	int fn_len = strlen(filename);
	char *ftype = (char *)malloc(5);
	for(int i = 0; i < 4; ++i) {
		ftype[i] = filename[fn_len - 4 + i];
	}

	ftype[4] = '\0';

	return ftype;
}

template <class int_type>
int_type bit_encode(const char c) {
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
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
    code_dict['C'] = bit_encode<int_type>('C');
    code_dict['G'] = bit_encode<int_type>('G');
    code_dict['T'] = bit_encode<int_type>('T');
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
void seq_decode(char* buf, const int len, const int_type seq_code, int_type* code_dict, const int_type b_mask) {
    for (int i=0;  i < len-1;  ++i) {
        const int_type b_code = (seq_code >> (bpb * (len - i - 2))) & b_mask;
        buf[i] = bit_decode<int_type>(b_code);
    }

    buf[len-1] = '\0';
}


template <class int_type>
void bit_load(const char* k_path, vector<char>& buffer, vector<tuple<int_type, int_type>>& k_vec, const int_type* code_dict, const int_type b_mask) {
    auto t_start = chrono_time();

    char* window = buffer.data();

    uintmax_t n_lines = 0;

    int fd;
    fd = open(k_path, O_RDONLY);

    int cur_pos = 0;
    int snp_pos = 0;

    char seq_buf[k];
	char snp_id[16];

	bool id_switch = false;
    bool has_wildcard = false;

    while (true) {

        const ssize_t bytes_read = read(fd, window, step_size);

        if (bytes_read == 0)
            break;

        if (bytes_read == (ssize_t) -1) {
            cerr << "unknown fetal error when reading " << k_path << endl;
            exit(EXIT_FAILURE);
        }

        for (int i = 0;  i < bytes_read;  ++i) {
            char c = toupper(window[i]);
            if (c == '\n') {
                ++n_lines;

                if (has_wildcard) {
                    has_wildcard = false;
                    continue;    
                }

                auto code = seq_encode<int_type>(seq_buf, k, code_dict, b_mask);

				snp_id[snp_pos] = '\0';
				int_type id_int = stoull(snp_id);

                k_vec.push_back(tuple<int_type, int_type>(code, id_int));

                cur_pos = 0;
				snp_pos = 0;

				id_switch = false;
            } else if (c == '\t'){
				id_switch = true;
			} else {
                if (c == 'N') {
                    has_wildcard = true;    
                }

				if (id_switch) {
					snp_id[snp_pos++] = c;	
				} else {
					seq_buf[cur_pos++] = c;
				}
            }
        }

    }

    auto timeit = chrono_time();
}

template <class int_type>
void binary_load(const char* k_path, vector<tuple<int_type, int_type>>& k_vec) {
	size_t filesize = get_fsize(k_path);
	// open file
	int fd = open(k_path, O_RDONLY, 0);
	assert(fd != -1);
	// execute mmap
	int_type* mmappedData = (int_type *) mmap(NULL, filesize, PROT_READ, MMAP_FLAGS, fd, 0);
	assert(mmappedData != MAP_FAILED);
	// write the mmapped data to stdout (= FD #1)
	// write(1, mmappedData, filesize);

	auto l_start = chrono_time();

	for (uint64_t i = 0; i < filesize/8; i=i+2) {
		auto kmer_int = mmappedData[i];
		auto snp = mmappedData[i+1];

		auto k_pair = make_tuple(kmer_int, snp);

		k_vec.push_back(k_pair);
	}

	//cleanup after mmap
	int rc = munmap(mmappedData, filesize);
	assert(rc == 0);
	close(fd);
}

template <class int_type>
bool cmp_triple_mid(const tuple<uint32_t, uint32_t, int_type> &a, const tuple<uint32_t, uint32_t, int_type> &b){
	    return get<1>(a) < get<1>(b);
}

template <class int_type>
bool cmp_triple_left(const tuple<uint32_t, uint32_t, int_type> &a, const tuple<uint32_t, uint32_t, int_type> &b){
	    return get<0>(a) < get<0>(b);
}

template <class int_type>
void split_binary_load(vector<tuple<int_type, int_type>>& k_vec, unordered_map<uint32_t, tuple<int64_t, int64_t>>& lmer_indx, vector<tuple<uint32_t, int_type>>& mmers, unordered_map<uint32_t, tuple<int64_t, int64_t>>& mmer_indx, vector<tuple<uint32_t, int_type>>& lmers) {
	auto l_start = chrono_time();

	cerr << chrono_time() << ":  " << "[OK] start to load DB " << endl;

	int64_t start = -1;
	int64_t end = -1;

	uint32_t cur_lmer = 2147483648;

	vector<tuple<uint32_t, uint32_t, int_type>> lm_vec;

	for (auto it = k_vec.begin(); it != k_vec.end(); ++it) {
		auto kmer_int = get<0>(*it);

		uint32_t lmer_int = (uint32_t)((kmer_int & 0xFFFFFFFF00000000LL) >> 32);
		uint32_t mmer_int = (uint32_t)(kmer_int & 0xFFFFFFFFLL);

		mmers.push_back(make_tuple(mmer_int, kmer_int));
		++end;

		auto lmer = lmer_int;

		auto lm_elem = make_tuple(lmer, mmer_int, kmer_int);
		lm_vec.push_back(lm_elem);

		if (cur_lmer != lmer) {
			if (start == -1){
				start = 0;
			} else {
				auto coord = make_tuple(start, end);
				lmer_indx.insert({cur_lmer, coord});

				start = end;
			}

			cur_lmer = lmer;
		}
	}

	if (lmer_indx.find(cur_lmer) == lmer_indx.end()) {
		auto coord = make_tuple(start, end);
		lmer_indx.insert({cur_lmer, coord});
	}

	sort(lm_vec.begin(), lm_vec.end(), cmp_triple_mid<int_type>);

	start = -1;
	end = -1;

	// using a very larger number as a starter
	auto cur_mmer = 2147483648;

	for (auto it = lm_vec.begin(); it != lm_vec.end(); ++it) {
		lmers.push_back(make_tuple(get<0>(*it), get<2>(*it)));
		auto mmer = get<1>(*it);

		++end;

		if (cur_mmer != mmer) {
            if (start == -1){
                start = 0;
            } else {
                auto coord = make_tuple(start, end);
                mmer_indx.insert({cur_mmer, coord});

                start = end;
            }

            cur_mmer = mmer;
        }
	}

	if (mmer_indx.find(cur_mmer) == mmer_indx.end()) {
		auto coord = make_tuple(start, end);
		mmer_indx.insert({cur_mmer, coord});
	}

	assert(lmers.size() == mmers.size());
	lm_vec.clear();
	cerr << chrono_time() << ":  " << "Done loading DB.  That took " << (chrono_time() - l_start) / 1000 << " seconds." << endl;
}


template <class int_type>
bool cmp_tuple(const tuple<int_type, int_type> &a, const tuple<int_type, int_type> &b){
	return get<0>(a) < get<0>(b);
}

template <class int_type>
void multi_kuniq(vector<string> kpaths, char* out_path, bool pool_sorted, bool suppress_output) {	
    int_type lsb = 1;
    int_type b_mask = (lsb << bpb) - lsb;

    int_type code_dict[1 << (sizeof(char) * 8)];
    make_code_dict<int_type>(code_dict);

	auto n_path = kpaths.size();

    vector<tuple<int_type, int_type>> kdb;
    vector<char> buffer(buffer_size);

	char* kdb_type = get_ftype(kpaths[0].c_str());

	if (strcmp(kdb_type, ".tsv") == 0) {
		bit_load<int_type>(kpaths[0].c_str(), buffer, kdb, code_dict, b_mask);
	} else if (strcmp(kdb_type, ".bin") == 0) {
		binary_load<int_type>(kpaths[0].c_str(), kdb);	
	} else {
		assert(false);
	}

    auto timeit = chrono_time();
	sort(kdb.begin(), kdb.end(), cmp_tuple<int_type>);
    cerr << "Sorting done! " << "It takes " << (chrono_time() - timeit) / 1000 << " secs" << endl;
    cerr << "the kmer list has " << kdb.size() << " kmers" << endl;
	cerr << "start to check conflicts" << endl;

    vector<tuple<int_type, int_type>> kpool;
    vector<tuple<int_type, int_type>> _kdb;

    for (int i = 1; i < n_path; ++i) {
		char* kp_type = get_ftype(kpaths[i].c_str());

		if (strcmp(kp_type, ".tsv") == 0) {
        	bit_load<int_type>(kpaths[i].c_str(), buffer, kpool, code_dict, b_mask);	
		} else if (strcmp(kp_type, ".bin") == 0) {
        	binary_load<int_type>(kpaths[i].c_str(), kpool);	
		} else {
			assert(false);
		}


		if (!pool_sorted) {
			sort(kpool.begin(), kpool.end(), cmp_tuple<int_type>);
		}

		auto ip = kpool.begin();
		for (auto it = kdb.begin(); it != kdb.end(); ++it) {
			while (get<0>(*it) > get<0>(*ip) && ip != kpool.end()) {
				++ip;
			}

			if(ip == kpool.end()) {
				_kdb.push_back(*it);
				continue;
			} 

			if (get<0>(*it) < get<0>(*ip)){
				_kdb.push_back(*it);	
			} 
		}

		kdb.swap(_kdb);
		cerr << "the kmer list has " << kdb.size() << " kmers after validation against: " << kpaths[i] << endl;

		_kdb.clear();
		kpool.clear();
    }


	vector<int_type> o_buff;

    for (auto it = kdb.begin(); it != kdb.end(); ++it) {
		o_buff.push_back(get<0>(*it));
		o_buff.push_back(get<1>(*it));
    }


	cerr << "the final kmer list has " << o_buff.size()/2<< " uniq kmers after purging conflicts" << endl;

	if (!suppress_output) {
		ofstream fh(out_path, ofstream::binary);

		fh.write((char*)&o_buff[0], o_buff.size() * sizeof(int_type));
		fh.close();
	}
}

void display_usage(char *fname){
	cout << "usage: " << fname << " -d db_path -o output_path [-S Sorted Pool Indicator] [-L pool_list_path] pool_path1 [pool_path2 ...] | \n";
}

int main(int argc, char** argv){		
	extern char *optarg; 
	extern int optind;

	bool dbflag = false;
	bool list_flag = false;
	bool sortflag = false;
	bool outsuppflag = false;

	char* fname = argv[0];
	char* db_path = (char *)"";
	char* list_path = (char *)"";
	char* opath = (char *)"/dev/stdout";


	int opt;
	while ((opt = getopt(argc, argv, "d:o:L:SOh")) != -1) {
		switch (opt) {
			case 'd':
				dbflag = true;
				db_path = optarg;
				break;
			case 'o':
				opath = optarg;
				break;
			case 'L':
				list_flag = true;
				list_path = optarg;
				break;
			case 'S':
				sortflag = true;
				break;
			case 'O':
				outsuppflag = true;
				break;
			case 'h': case '?':
				display_usage(fname);
				exit(1);
		}	
	}

	cout << fname << '\t' << db_path << '\t' << endl;

	if (!dbflag) {
		cout << "missing argument: -d <sckmerdb_path: string>\n";
		display_usage(fname);
		exit(1);
	}

	if (list_flag) {
		cout << "program reads a list of kmer pools for checking kmer uniqueness: " << list_path << endl;
	} else {
		int in_pos = optind;
		if (optind == argc) {
			cout << "missing argument: input (>1)\n";
			display_usage(fname);
			exit(1);
		}
	}

	vector<string> kp_paths;
	
	kp_paths.push_back(string(db_path));
	
	if (list_flag) {
		ifstream file(list_path);
		string line;
		while (getline(file, line)) {
			string tmp_line = line;
			kp_paths.push_back(tmp_line);
		}
	}

	for(; optind < argc; optind++) {
		string kp_path(argv[optind]);
		kp_paths.push_back(kp_path);
	}

	cout << "the number of kmer pool files that will be used: " << kp_paths.size() - 1 << endl;

	if (argc == 2 && string(argv[1]) == "-h") {
		display_usage(argv[0]);
	} else if (argc >= 2) {
		multi_kuniq<uint64_t>(kp_paths, opath, sortflag, outsuppflag);		
	} else {
        cerr << argv[0] << " takes at least two arguments!" << endl;
		display_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    return 0;
}

