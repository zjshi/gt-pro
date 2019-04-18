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

template <class int_type>
void bit_load(const char* k_path, vector<tuple<int_type, int_type>>& k_vec, const int_type* code_dict, const int_type b_mask) {

    auto t_start = chrono_time();
    assert(k_path);
    cerr << "Loading file " << k_path << "." << endl;

    FILE *fp = fopen(k_path, "r");
    if (fp == NULL) {
        cerr << "Trouble opening file " << k_path << "." << endl;
        exit(EXIT_FAILURE);
    }

    // 3 columns in each line.
    struct Column {
        char* str;
        size_t size;
        ssize_t len;
        int sep;
        Column(int separator) : str(NULL), size(0), len(0), sep(separator) {};
        ~Column() {
            if (str) {
                free(str);
            }
        }
        inline void read_from_file(FILE* fp) {
            // getdelim will (re)allocate memory for str as needed.
            len = getdelim(&str, &size, sep, fp);
            if (len > 0 && str[len - 1] == sep) {
                str[--len] = '\0';
            }
        }
        inline bool is_numeric() {
            char* s = str;
            while (isdigit(*s)) {
                ++s;
            }
            return *s == '\0';
        }
        inline bool does_not_contain_wildcard() {
            char* s = str;
            while (*s && *s != 'n' && *s != 'N') {
                ++s;
            }
            return *s == '\0';
        }
    };

    Column cols[] = { Column('\t'), Column('\t'), Column('\n') };

    while (true) {
        for (auto& c : cols) {
            c.read_from_file(fp);
        }
        const bool end_of_file = cols[0].len < 0;
        const bool some_columns_empty = cols[0].len <= 0 || cols[1].len <= 0 || cols[2].len <= 0;
        assert(some_columns_empty == end_of_file && "Every line should have 3 non-empty columns, separated by tabs.");
        if (end_of_file) {
            break;
        }
        assert(cols[0].len == K && "First column needs to be a K-mer (contain exactly K characters)");
        assert(cols[1].len <= 2 && cols[1].is_numeric() && "Second column needs to be decimal with at most 2 digits.");
        assert(cols[2].len <= 16 && cols[2].is_numeric() && "Third and last column needs to be decimal with at most 16 digits.");
        auto snp_offset = strtoul(cols[1].str, NULL, 10);
        assert(snp_offset >= 0 && snp_offset < K && "Second column needs to be decimal in range 0 .. K-1, inclusive.");
        int_type snp = strtoull(cols[2].str, NULL, 10);
        if (cols[0].does_not_contain_wildcard()) {
            auto kmer = seq_encode<int_type, K>(cols[0].str, code_dict, b_mask);
            k_vec.push_back(tuple<int_type, int_type>(kmer, snp));
        }
    }
    // Done parsing file.
    cerr << "Loaded file " << k_path << " in " << (chrono_time() - t_start) / 1000.0 << " secs." << endl;
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
    // We want the most signifficant 6 digits.  The number of digits to begin with is at most 16.
    // Here is a slow way to do it:
    //     x = stoi(to_string(get_snp(kmer_snp)).substr(0, 6));
    // The fast way is below.  This actually makes a big difference in overall perf.
    auto x = get_snp(kmer_snp);
    // if it has 13 or more digits, remove the last 7
    if (x >= 1000000000000) { // 10**12
        x /= 10000000;        // 10**7
    }
    // if it has 10 or more digits, remove the last 4
    if (x >= 1000000000) {   // 10**9
        x /= 10000;          // 10**4
    }
    // so long as it still has 7 or more digits, keep chopping off the last digit
    while (x >= 1000000) {   // 10**6
        x /= 10;             // 10**1
    }
    return x;
}

template <class int_type>
void multi_btc64(int n_path, char** kpaths) {
    int_type lsb = 1;
    int_type b_mask = (lsb << BITS_PER_BASE) - lsb;

    auto code_dict = new_code_dict<int_type>();

    vector<tuple<int_type, int_type>> kdb;

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
        // cerr << "consider kmer " << hex << kmer << " species " << dec << species << endl;
        bool kmer_is_monospecific = true;
        auto next = current + 1;
        while (next != kdb.end() && get_kmer(*next) == kmer) {
            // cerr << "consider next species " << get_species(*next);
            kmer_is_monospecific = kmer_is_monospecific && (get_species(*next) == species);
            // cerr << " monospecific " << kmer_is_monospecific << endl;
            ++next;
        }
        if (kmer_is_monospecific) {
            // Every snp for this kmer is from the same species.  Emit.
            for (auto kmer_snp = current;  kmer_snp != next;  ++kmer_snp) {
                // assert(get_kmer(*kmer_snp) == kmer);
                o_buff.push_back(kmer);
                o_buff.push_back(get_snp(*kmer_snp));
                // cerr << "emit " << hex << kmer << " " << dec << get_snp(*kmer_snp) << endl;
            }
            // Count stats.
            monospecies_kmers += (next - current);
            ++monospecies_unique_kmers;
        } else {
            // The current kmer hits multiple species.  Suppress.
            // Just count some stats.
            multispecies_kmers += (next - current);
            ++multispecies_unique_kmers;
            // cerr << "do not emit " << hex << kmer << endl;
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
