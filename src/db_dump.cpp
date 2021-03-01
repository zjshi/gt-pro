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


// this program scans its input (fastq text stream) for forward k mers,

// usage:
//    g++ -O3 --std=c++11 -o vfkmrz_bunion vfkmrz_bunion.cpp
//    ./vfkmrz_bunion -k1 </path/to/kmer_list1> -k2 </path/to/kmer_list2>
//
// standard fastq format only for input, otherwise failure is almost guaranteed. 

// global variable declaration starts here
constexpr auto k = 31;

// set operation mode
// valid values: 0, 1, 2
// 0 is set union operation; 1 is set intersection operation; 2 is set difference([set1-set2]);
constexpr auto s_mod = 0;

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
void seq_decode(char* buf, const int len, const int_type seq_code, const int_type* code_dict, const int_type b_mask) {
    for (int i=0;  i < len-1;  ++i) {
        const int_type b_code = (seq_code >> (bpb * (len - i - 2))) & b_mask;
        buf[i] = bit_decode<int_type>(b_code);
    }

    buf[len-1] = '\0';
}


template <class int_type>
void bit_load(const char* k_path, vector<int_type>& buffer, vector<int_type>& k_vec, const int_type* code_dict, const int_type b_mask) {
    auto t_start = chrono_time();

    int_type* window = buffer.data();

    uintmax_t n_lines = 0;

    int cur_pos = 0;
    char seq_buf[k+1];

    //auto fh = fstream(out_path, ios::out | ios::binary);

    int fd;
    fd = open(k_path, O_RDONLY);

    bool has_wildcard = false;

    ofstream fh(out_path, ofstream::binary);

    while (true) {

        const ssize_t bytes_read = read(fd, window, step_size);

        if (bytes_read == 0)
            break;

        if (bytes_read == (ssize_t) -1) {
            cerr << "unknown fetal error!" << endl;
            exit(EXIT_FAILURE);
        }

        for (int i = 0;  i < bytes_read / 8;  ++i) {
            int_type code = window[i];
                // auto code = seq_encode<int_type>(seq_buf, k, code_dict, b_mask);

			if ( i % 2 == 0) {
				seq_decode(seq_buf, k+1, code, code_dict, b_mask);    
				fh << seq_buf << '\t';
			} else {
				fh << code << '\n';
			}
        }

    }

    close(fd);
    fh.close();

    auto timeit = chrono_time();
}

template <class int_type>
void multi_btdc64(int n_path, char** kpaths) {	
    int_type lsb = 1;
    int_type b_mask = (lsb << bpb) - lsb;

    int_type code_dict[1 << (sizeof(char) * 8)];
    make_code_dict<int_type>(code_dict);
	
	/*
	uint32_t lsb_32 = 1;
    uint32_t b_mask_32 = (lsb_32 << bpb) - lsb_32;

    uint32_t code_dict_32[1 << (sizeof(char) * 4)];
    make_code_dict<uint32_t>(code_dict_32);
	*/


    vector<int_type> kdb;
    vector<int_type> buffer(buffer_size);

    for (int i = 1; i < n_path; ++i) {
        cerr << kpaths[i] << endl;
        bit_load<int_type>(kpaths[i], buffer, kdb, code_dict, b_mask);	
    }

    auto timeit = chrono_time();
}

int main(int argc, char** argv){		
    if (argc >= 2) {
        multi_btdc64<uint64_t>(argc, argv);		
    } else {
        cerr << argv[0] << "takes at least one arguments!" << endl;
        exit(EXIT_FAILURE);
    }

    return 0;
}
