#ifndef INDEL_DETECTION_RM_INDEX_H
#define INDEL_DETECTION_RM_INDEX_H

#include "../common/packed_db.h"

#define MAX_KMER_CNT 512 //1000

class ReferenceMappingLookupTable
{
public:
    ReferenceMappingLookupTable() : kmer_starts(NULL), kmer_cnts(NULL), offset_list(NULL) {}
    ~ReferenceMappingLookupTable()
    {
        if (kmer_starts) delete[] kmer_starts;
        if (kmer_cnts) delete[] kmer_cnts;
        if (offset_list) delete[] offset_list;
    }

    void build(PackedDB* ref, const int kmer_size, const int num_threads, const idx KmerCntCutoff);

    bool extract_kmer_list(const u32 hash, i64*& offset_arr, int& arr_size)
    {
        offset_arr = offset_list + kmer_starts[hash];
        arr_size = kmer_cnts[hash];
        return arr_size;
    }

private:
    

private:
    i64*   kmer_starts;
    short*  kmer_cnts;
    i64*   offset_list;
};

typedef ReferenceMappingLookupTable LookupTable;

#endif //INDEL_DETECTION_RM_INDEX_H
