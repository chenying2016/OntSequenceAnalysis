#ifndef LOOKUP_TABLE_H
#define LOOKUP_TABLE_H

#include "packed_db.h"

class LookupTable
{
public:
    typedef u64 index_t;

public:
    LookupTable(): kmer_cnts(NULL), kmer_starts(NULL), offset_list(NULL) {}

    ~LookupTable() {
        if (kmer_cnts) delete[] kmer_cnts;
        if (kmer_starts) delete[] kmer_starts;
        if (offset_list) delete[] offset_list;
    }

    void build(const PackedDB* ref, const int kmer_size, const int kmer_cnt_cutoff, const int num_threads);

    bool extract_kmer_list(const u32 hash, index_t*& offsets, int& n) {
        offsets = kmer_starts[hash];
        n = kmer_cnts[hash];
        return n;
    }

private:
    short*      kmer_cnts;
    index_t**   kmer_starts;
    index_t*    offset_list;
};

#endif // LOOKUP_TABLE_H
