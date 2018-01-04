#include "lookup_table.h"
#include "lookup_table_bucket_sort.h"
    

void 
LookupTable::build(const PackedDB* ref, const int kmer_size, const int kmer_cnt_cutoff, const int num_threads)
{
    ns_mc_lookup_table::build_lookup_table(ref, 
            kmer_size,
            kmer_cnt_cutoff,
            num_threads,
            kmer_cnts,
            kmer_starts,
            offset_list);
}
