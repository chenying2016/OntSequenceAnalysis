#include "lookup_table.h"
#include "timer.h"

#include <algorithm>

using namespace std;

//static const int KmerCntCutoff = MAX_KMER_CNT;

#define get_next_hash { \
	const u8 c = ref->get_char(offset + j); \
	hash = ((hash << 2) | c) & hash_mask; \
}
	

void 
ReferenceMappingLookupTable::build(PackedDB* ref, const int kmer_size, const int num_threads , const idx KmerCntCutoff)
{
	DynamicTimer dt(__func__);
	
	u64 max_hash = 1; max_hash = max_hash << (kmer_size * 2);
	u64 hash_mask = max_hash - 1;
	kmer_cnts = new short[max_hash];
	fill(kmer_cnts, kmer_cnts + max_hash, 0);
	int* cnt_table = new int[KmerCntCutoff + 2];
	for (int i = 0; i <= KmerCntCutoff; ++i) cnt_table[i] = i + 1;
	cnt_table[KmerCntCutoff + 1] = KmerCntCutoff + 1;
	
	const idx ns = ref->num_seqs();
	for (idx i = 0; i < ns; ++i) {
		const idx offset = ref->seq_offset(i);
		const idx size = ref->seq_size(i);
		u64 hash = 0;
		for (idx j = 0; j < kmer_size - 1; ++j) {
			get_next_hash;
		}
		for (idx j = kmer_size - 1; j < size; ++j) {
			get_next_hash;
			kmer_cnts[hash] = cnt_table[kmer_cnts[hash]];
		}
	}
	
	idx num_kmers = 0;
	kmer_starts = new i64[max_hash];
	for (u64 i = 0; i != max_hash; ++i) {
		if (kmer_cnts[i] <= KmerCntCutoff) {
			kmer_starts[i] = num_kmers;
			num_kmers += kmer_cnts[i];
		} else {
			kmer_starts[i] = -1;
		}
		kmer_cnts[i] = 0;
	}
	
	offset_list = new i64[num_kmers];
	i64 num_added = 0;
	for (idx i = 0; i < ns; ++i) {
		const idx offset = ref->seq_offset(i);
		const idx size = ref->seq_size(i);
		u64 hash = 0;
		for (idx j = 0; j < kmer_size - 1; ++j) {
			get_next_hash;
		}
		for (idx j = kmer_size - 1; j < size; ++j) {
			get_next_hash;
			if (kmer_starts[hash] != -1) {
				idx k = kmer_starts[hash] + kmer_cnts[hash];
				offset_list[k] = offset + j + 1 - kmer_size;
				++kmer_cnts[hash];
				++num_added;
			} else {
				r_assert(kmer_cnts[hash] == 0);
			}
		}
	}
	r_assert(num_added == num_kmers);
	delete[] cnt_table;
}
