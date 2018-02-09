#ifndef WORD_FINDER_H
#define WORD_FINDER_H

#include "../klib/kvec.h"
#include "../common/map_options.h"
#include "../common/ontcns_aux.h"
#include "../common/packed_db.h"
#include "../lookup_table/lookup_table.h"
#include "chain_dp.h"
#include "word_finder_aux.h"

typedef struct
{
	vec_u64				hash_list;
	vec_chain_seed		seed_list;
	vec_chain_seed		chain_seed_list;
	vec_int				ddf_scores;
	vec_can				candidates;
	vec_bseeds_index	block_seeds_index;
	ChainDpData*		chain_data;
} WordFindData;

WordFindData*
new_WordFindData(int kmer_size, int block_score_cutoff);

WordFindData*
free_WordFindData(WordFindData* data);

void
find_candidates(const char* read,
				const int read_size,
				const int qid,
				const int qdir,
				const int read_start_id,
				const int reference_start_id,
				BOOL pairwise,
				PackedDB* reference,
				LookupTable* lktbl,
				MapOptions* options,
				WordFindData* wfdata,
			    vec_can* candidates);

#endif // WORD_FINDER_H
