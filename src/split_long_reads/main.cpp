#include "../common/fasta_reader.h"

#include <fstream>

using namespace std;

void
print_usage(const char* prog) 
{
    const char sep = ' ';
    cerr << "USAGE:" << endl
         << prog << sep
         << "min-length" << sep
         << "max-length" << sep
         << "input" << sep
         << "output" << endl;
}

void
split_read(Sequence& read, const idx min_size, const idx max_size, idx& read_id, ostream& out)
{
    idx s = read.size();
    if (s < min_size) return;
    if (max_size < 0 || s < max_size) {
        out << '>' << "ontcns_" << read_id << '\n';
        read_id++;
        const Sequence::str_t& seq = read.sequence();
        for (size_t i = 0; i != seq.size(); ++i) {
            out << seq[i];
        }
        out << '\n';
        return;
    }
    
    idx L = 0, R = 0;
    const Sequence::str_t& seq = read.sequence();
    while (L < s) {
        R = min(L + max_size, s);
        if (R - L < min_size) break;
        out << '>' << "ontcns_" << read_id << '\n';
        ++read_id;
        for (idx i = L; i < R; ++i) {
            out << seq[i];
        }
        out << '\n';
        L = R;
    }
}

int main(int argc, char* argv[])
{
    if (argc != 5) {
        print_usage(argv[0]);
        return 1;
    }
    idx min_size = atoll(argv[1]);
    idx max_size = atoll(argv[2]);
    const char* input = argv[3];
    const char* output = argv[4];
    Sequence read;
    FastaReader reader(input);
    idx read_id = 0;
    ofstream out(output);

    while (1) {
        idx s = reader.read_one_seq(read);
        if (s == -1) break;
        split_read(read, min_size, max_size, read_id, out);
    }

    out.close();
    return 0;
}
