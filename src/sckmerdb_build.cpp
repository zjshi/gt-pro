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

// global variable declaration starts here
constexpr auto K = 31;

// parameters for <unistd.h> file read; from the source of GNU coreutils wc
constexpr auto STEP_SIZE = 256 * 1024 * 1024;
constexpr auto BUFFER_SIZE = 256 * 1024 * 1024;

// output file path
constexpr auto OUT_PATH = "/dev/stdout";

// get time elapsed since when it all began in milliseconds.
long chrono_time() {
    using namespace chrono;
    return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

constexpr int BITS_PER_BASE = 2;

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
vector<int_type> new_code_dict() {

    constexpr auto CHAR_LIMIT = 1 << (sizeof(char) * 8);
    vector<int_type> code_dict;
    for (uint64_t c = 0; c < CHAR_LIMIT; ++c) {
        // This helps us detect non-nucleotide characters on encoding.
        code_dict.push_back(-1);
    }

    code_dict['A'] = code_dict['a'] = bit_encode<int_type>('A');
    code_dict['C'] = code_dict['c'] = bit_encode<int_type>('C');
    code_dict['G'] = code_dict['g'] = bit_encode<int_type>('G');
    code_dict['T'] = code_dict['t'] = bit_encode<int_type>('T');

    return code_dict;
}

template <class int_type, int len>
int_type seq_encode(const char* buf, const int_type* code_dict, const int_type b_mask) {
    int_type seq_code = 0;
    for (int i=0;  i < len;  ++i) {
        const int_type b_code = code_dict[buf[i]];
        assert(b_code != -1 && "Trying to encode a character that is not ACTG or actg.");
        seq_code |= ((b_code & b_mask) << (BITS_PER_BASE * (len - i - 1)));
    }
    return seq_code;
}

bool is_numeric(char* str) {
    while (isdigit(*str)) {
        ++str;
    }
    return *str == '\0';
}

template <class int_type>
void bit_load(const char* k_path, vector<tuple<int_type, int_type>>& k_vec, const int_type* code_dict, const int_type b_mask) {

    assert(k_path);
    cerr << "Loading file " << k_path << "." << endl;

    FILE *fp = fopen(k_path, "r");
    if (fp == NULL) {
        cerr << "Trouble opening file " << k_path << "." << endl;
        exit(EXIT_FAILURE);
    }

    // The getline function (re)allocates memory for the line as needed.
    char *line = NULL;
    size_t buf_size = 0;
    ssize_t line_len;
    auto t_start = chrono_time();

    while ((line_len = getline(&line, &buf_size, fp)) != -1) {
        // A line of line_len characters has been read.   Squash the terminating \n, if present.
        if (line_len >= 1 && line[line_len - 1] == '\n') {
            line[line_len - 1] = '\0';
        }
        char *first_tab;
        char *second_tab;
        int_type kmer;
        unsigned int snp_offset;
        int_type snp;
        bool has_wildcard;
        // Parse first column: the K-mer.
        {
            // Identify extents of first_column, and make first_column into a 0-terminated string.
            char *first_column = line;  // note: getline may move line to a different address in each iteration of this loop
            first_tab = index(first_column, '\t');
            assert(first_tab && ((first_tab - first_column) == K) && "First column needs to be a K-mer (contain exactly K characters).");
            *first_tab = '\0';  // now first_column is a 0-terminated string.
            has_wildcard = index(first_column, 'N') || index(first_column, 'n');
            // This now supports both lowercase and uppercase ACTG characters.
            if (!(has_wildcard)) {
                kmer = seq_encode<int_type, K>(first_column, code_dict, b_mask);
            }
        }
        // Parse second column: the snp offset within the K-mer (integer range 0..K-1).
        {
            char * second_column = first_tab + 1;
            second_tab = index(second_column, '\t');
            assert(second_tab && "Line must contain exactly 3 columns.");
            *second_tab = '\0';  // now second_column is a 0-terminated string.
            assert(is_numeric(second_column) && "Second column needs to be decimal in range 0 .. K-1, inclusive.");
            assert(strlen(second_column) <= 5 && "Second column needs to be decimal in range 0 .. K-1, inclusive.");
            snp_offset = strtoul(second_column, NULL, 10);
            assert(snp_offset >= 0 && "Second column needs to be decimal in range 0 .. K-1, inclusive.");
            assert(snp_offset < K && "Second column needs to be decimal in range 0 .. K-1, inclusive.");
        }
        // Parse third and last column: decimal literal encoding species, major/minor allele, and genomic position of SNP.
        {
            char* third_column = second_tab + 1;
            assert(is_numeric(third_column) && "Third and last column needs to be decimal literal with at most 16 digits.");
            assert(strlen(third_column) <= 16 && "Third and last column needs to be decimal literal with at most 16 digits.");
            snp = strtoull(third_column, NULL, 10);
        }
        // Finally, emit the (kmer, snp) pair.  The snp_offset is still ignored at the moment.
        if (!(has_wildcard)) {
            k_vec.push_back(tuple<int_type, int_type>(kmer, snp));
        }
    }
    // Done parsing file.
    free(line);
    cerr << "Loaded file " << k_path << " in " << (chrono_time() - t_start) / 1000 << " secs." << endl;
}

template <class int_type>
int_type get_kmer(const tuple<int_type, int_type>& kmer_snp) {
    return get<0>(kmer_snp);
}

template <class int_type>
int_type get_snp(const tuple<int_type, int_type>& kmer_snp) {
    return get<1>(kmer_snp);
}

template <class int_type>
int_type get_species(const tuple<int_type, int_type>& kmer_snp) {
    return stoi(to_string(get_snp(kmer_snp)).substr(0, 6));
}

template <class int_type>
void multi_btc64(int n_path, char** kpaths) {
    int_type lsb = 1;
    int_type b_mask = (lsb << BITS_PER_BASE) - lsb;

    auto code_dict = new_code_dict<int_type>();

    vector<tuple<int_type, int_type>> kdb;
    vector<char> buffer(BUFFER_SIZE);

    if (n_path == 0) {
        cerr << "No paths specified on command line --> will read from /dev/stdin." << endl;
        n_path = 1;
        char* stdin = "/dev/stdin";
        kpaths = &stdin;
    }

    auto timeit = chrono_time();

    for (int i = 0; i < n_path; ++i) {
        bit_load<int_type>(kpaths[i], kdb, code_dict.data(), b_mask);
    }
    cerr << "Loaded all files! " << "It took " << (chrono_time() - timeit) / 1000 << " secs." << endl;
    cerr << "The unfiltered kmer list has " << kdb.size() << " kmers" << endl;

    timeit = chrono_time();

    sort(kdb.begin(), kdb.end());
    cerr << "Sorting done! " << "It took " << (chrono_time() - timeit) / 1000 << " secs." << endl;

    char seq_buf[K + 1];

    uint64_t multispecies_kmers = 0;
    uint64_t multispecies_unique_kmers = 0;
    uint64_t monospecies_kmers = 0;
    uint64_t monospecies_unique_kmers = 0;

    vector<int_type> o_buff;

    cerr << "Starting to check conflicts:  Will filter out kmers that occur in multiple species." << endl;

    timeit = chrono_time();

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
            for (auto kmer_snp = current;  kmer_snp != next;  ++kmer_snp) {
                assert(get_kmer(*kmer_snp) == kmer);
                o_buff.push_back(kmer);
                o_buff.push_back(get_snp(*kmer_snp));
            }
            // Count stats.
            monospecies_kmers += (next - current);
            ++monospecies_unique_kmers;
        } else {
            // The current kmer hits multiple species.  Suppress.
            // Just count some stats.
            multispecies_kmers += (next - current);
            ++multispecies_unique_kmers;
        }
        current = next;
    }

    cerr << "Filtering done! " << "It took " << (chrono_time() - timeit) / 1000 << " secs." << endl;
    cerr << "The filtered kmer list has " << o_buff.size()/2<< " kmers after purging conflicts." << endl;
    cerr << "Monospecies kmers: " <<  monospecies_kmers << " (" << monospecies_unique_kmers << " unique)." << endl;
    cerr << "Multispecies kmers: " <<  multispecies_kmers << " (" << multispecies_unique_kmers << " unique)." << endl;

    ofstream fh(OUT_PATH, ofstream::binary);

    fh.write((char*)&o_buff[0], o_buff.size() * sizeof(int_type));
    fh.close();
}

void display_usage(char *fname){
    cout << "usage: " << fname << " fpath [fpath ...]\n";
}

int main(int argc, char** argv){
    if (argc == 2 && (string(argv[1]) == "-h")) {
        display_usage(argv[0]);
    } else {
        multi_btc64<uint64_t>(argc - 1, argv + 1);
    }
    return 0;
}
