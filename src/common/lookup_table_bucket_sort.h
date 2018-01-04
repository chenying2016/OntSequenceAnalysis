#ifndef LOOKUP_TABLE_BUCKET_SORT_H
#define LOOKUP_TABLE_BUCKET_SORT_H

#include "packed_db.h"

namespace ns_mc_lookup_table {

void 
build_lookup_table(const PackedDB* ref, 
				   const int kmer_size,
				   const int kmer_cnt_cutoff,
				   const int num_threads,
				   short*& kmer_cnts, 
				   u64**& kmer_starts, 
				   u64*& offset_list);

} // ns_mc_lookup_table

#endif // LOOKUP_TABLE_BUCKET_SORT_H
