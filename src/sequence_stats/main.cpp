#include "../common/fasta_reader.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <numeric>

using namespace std;

struct LengthStats
{
    idx num_sequences;
    idx max_length;
    idx min_length;
    idx total_length;
    idx avg_length;
    idx median_length;
    idx n25_length;
    idx n25_cnt;
    idx n50_length;
    idx n50_cnt;
    idx n75_length;
    idx n75_cnt;
};

void print_length(idx length) {
    ostringstream os;
    os << length;
    string str = os.str();
    size_t n = str.size();
    int c = 0;
    char buf[256];
    char* p = buf;
    for (size_t i = n; i > 0; --i) {
        *p++ = str[n - 1];
        ++c;
        if (c == 3) {
            *p++ = ',';
            c = 0;
        }
    }
    while (p > buf) {
        cout << *(p - 1);
        --p;
    }
    cout << endl;
}

struct IdxGT
{
    bool operator()(const idx a, const idx b) {
        return a > b;
    }
};

void
read_sequence_length(const char* path, vector<idx>& length, const idx length_cutoff)
{
    FastaReader reader(path);
    Sequence seq;
    while (1) {
        idx s = reader.read_one_seq(seq);
        if (s == -1) break;
		if (s < length_cutoff) continue;
        length.push_back(s);
    }
    sort(length.begin(), length.end(), IdxGT());
}

void
calc_nxx_stats(vector<idx>& length, const idx nxx_sum, idx& nxx_length, idx& nxx_cnt)
{
    nxx_cnt = 0;
    idx sum = 0;
    for (size_t i = 0; i != length.size(); ++i) {
        ++nxx_cnt;
        sum += length[i];
        if (sum >= nxx_sum) {
            nxx_length = length[i];
            break;
        }
    }
}

void 
calc_length_stats(vector<idx>& length, LengthStats& stats)
{
    if (length.empty()) return;
    size_t n = length.size();
    stats.num_sequences = static_cast<idx>(length.size());
    stats.max_length = length[0];
    stats.min_length = length[n - 1];
    stats.median_length = length[n / 2];

    idx sum = 0;
    for (size_t i = 0; i < n; ++i) sum += length[i];
    idx n25_sum = sum / 4;
    idx n50_sum = n25_sum * 2;
    idx n75_sum = n25_sum * 3;
    stats.total_length = sum;
    stats.avg_length = sum / stats.num_sequences;
    calc_nxx_stats(length, n25_sum, stats.n25_length, stats.n25_cnt);
    calc_nxx_stats(length, n50_sum, stats.n50_length, stats.n50_cnt);
    calc_nxx_stats(length, n75_sum, stats.n75_length, stats.n75_cnt);
}

void 
print_length_stats(const LengthStats& stats)
{
    cerr << "Number of sequences:\t\t" << stats.num_sequences << endl;
    cerr << "Number of bps:\t\t\t" << stats.total_length << endl;
    cerr << "Max length:\t\t\t" << stats.max_length << endl;
    cerr << "Min length:\t\t\t" << stats.min_length << endl;
    cerr << "Average length:\t\t\t" << stats.avg_length << endl;
    cerr << "Median length:\t\t\t" << stats.median_length << endl;
    cerr << "N25 stats:\t\t\t" << "25% of total sequence is contained in the " << stats.n25_cnt << " sequences >= " << stats.n25_length << " bp" << endl;
    cerr << "N50 stats:\t\t\t" << "50% of total sequence is contained in the " << stats.n50_cnt << " sequences >= " << stats.n50_length << " bp" << endl;
    cerr << "N50 stats:\t\t\t" << "75% of total sequence is contained in the " << stats.n75_cnt << " sequences >= " << stats.n75_length << " bp" << endl;
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
        cerr << "USAGE:" << endl
             << argv[0] << ' ' << "sequences-name" << ' ' << "sequence-length-cutoff" << endl;
        return 1;
    }

    vector<idx> length;
    LengthStats stats;
	const char* sequence_path = argv[1];
	const idx sequence_length_cutoff = atoll(argv[2]);
    read_sequence_length(sequence_path, length, sequence_length_cutoff);
    calc_length_stats(length, stats);;
    print_length_stats(stats);
    return 0;
}
