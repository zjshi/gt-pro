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


// this program simply scans its input (text-based canonical k-mers) and 
// convert them into a binary format compatible with db_uniq 
// and print output to stdout

// usage:
//    g++ -O3 --std=c++11 -o mk_pool mk_pool.cpp
//    ./mk_pool <binary input>
//

// global variable declaration starts here
constexpr auto k = 31;

// parameters for <unistd.h> file read; from the source of GNU coreutils wc
constexpr auto step_size = 256 * 1024 * 1024;
constexpr auto buffer_size = 256 * 1024 * 1024;

// output file path
constexpr auto out_path = "/dev/stdout";

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

void make_comp_map(char* comp_map) {
    comp_map['A'] = 'T';
    comp_map['C'] = 'G';
    comp_map['G'] = 'C';
    comp_map['T'] = 'A';
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

// this function reads KMC k-mers, 
// then generates reverse complementary k-mers,
// and encode them into 64-bit integers
// and store them in a binary form:
// interger1[64-bit; k-mer 1],interger2[64-bit; snp ID 1],interger3[64-bit; k-mer 2],interger4[64-bit; snp ID 2]...
template <class int_type>
void bit_load(const char* k_path, vector<char>& buffer, vector<tuple<int_type, int_type>>& k_vec, const int_type* code_dict, const int_type b_mask) {
    auto t_start = chrono_time();

	char comp_map[1 << (sizeof(char) * 8)];
    make_comp_map(comp_map);

    char* window = buffer.data();

    uintmax_t n_lines = 0;

    int fd;
    fd = open(k_path, O_RDONLY);

    int cur_pos = 0;
    int snp_pos = 0;

    char seq_buf[k];
    char rcseq_buf[k];
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

				for (int l = k-1; l >= 0; --l) {
					rcseq_buf[k-1-l] = comp_map[seq_buf[l]];
				}

                auto code1 = seq_encode<int_type>(seq_buf, k, code_dict, b_mask);
				auto code2 = seq_encode<int_type>(rcseq_buf, k, code_dict, b_mask);

				snp_id[snp_pos] = '\0';
				int_type id_int = stoull(snp_id);

                k_vec.push_back(tuple<int_type, int_type>(code1, id_int));
                k_vec.push_back(tuple<int_type, int_type>(code2, id_int));

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
bool cmp_tuple(const tuple<int_type, int_type> &a, const tuple<int_type, int_type> &b){
	return get<0>(a) < get<0>(b);
}

template <class int_type>
void multi_btc64(int n_path, char** kpaths) {	
    int_type lsb = 1;
    int_type b_mask = (lsb << bpb) - lsb;

    int_type code_dict[1 << (sizeof(char) * 8)];
    make_code_dict<int_type>(code_dict);

    vector<tuple<int_type, int_type>> kdb;
    vector<char> buffer(buffer_size);

    for (int i = 1; i < n_path; ++i) {
        cerr << kpaths[i] << endl;

		char* kp_type = get_ftype(kpaths[i]);
		
		bit_load<int_type>(kpaths[i], buffer, kdb, code_dict, b_mask);	
    }

    auto timeit = chrono_time();

	sort(kdb.begin(), kdb.end(), cmp_tuple<int_type>);
    cerr << "Sorting done! " << "It takes " << (chrono_time() - timeit) / 1000 << " secs" << endl;
    cerr << "the kmer list has " << kdb.size() << " kmers" << endl;

	vector<int_type> o_buff;
	bool eq_last = false;

    cerr << "start to check conflicts" << endl;
    for (auto it = kdb.begin(); it+1 != kdb.end(); ++it) {
		if (get<0>(*it) == get<0>(*(it+1))) {
			eq_last = true;
		} else {
			if (eq_last) {
				eq_last = false;
			} else {
				o_buff.push_back(get<0>(*it));
				o_buff.push_back(get<1>(*it));
			}
		}
    }

	if (!eq_last) {
		auto last_elem = kdb.back();
		o_buff.push_back(get<0>(last_elem));
		o_buff.push_back(get<1>(last_elem));
	}

    ofstream fh(out_path, ofstream::binary);

    fh.write((char*)&o_buff[0], o_buff.size() * sizeof(int_type));
    fh.close();
}

void display_usage(char *fname){
	cout << "usage: " << fname << " fpath [fpath ...]\n";
}

int main(int argc, char** argv){		
	if (argc == 2 && string(argv[1]) == "-h") {
		display_usage(argv[0]);
	} else if (argc >= 2) {
        multi_btc64<uint64_t>(argc, argv);		
    } else {
        cerr << argv[0] << " reads from stdin or takes at least one arguments!" << endl;
		display_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    return 0;
}
