#include "word_finder.h"

#include <math.h>

#include "../common/ksort.h"
#include "../common/kstring.h"
#include "../common/kvec.h"
#include "../common/packed_db.h"
#include "../lookup_table/lookup_table.h"
#include "word_finder_aux.h"

KSORT_INIT(ChainSeed, ChainSeed, ChainSeedLT)

double
ddf_score(const size_t qoff1, const size_t soff1, const size_t qoff2, const size_t soff2)
{
	if (qoff1 == qoff2 || soff1 == soff2) return 1.0;
	double qdiff = 1.0 * ((qoff1 > qoff2) ? (qoff1 - qoff2) : (qoff2 - qoff1));
	double sdiff = 1.0 * ((soff1 > soff2) ? (soff1 - soff2) : (soff2 - soff1));
	double s = sdiff / qdiff;
	s = fabs(1.0 - s);
	return s;
}

int
extract_hash_values(const char* read,
        const int read_size,
        const int kmer_size,
        const int scan_window,
        vec_u64* hash_list)
{
    kv_clear(*hash_list);
    for (int i = 0; i <= read_size - kmer_size; i += scan_window) {
        u64 hash = 0;
        for (int j = 0; j < kmer_size; ++j) {
            const u8 c = (u8)read[i + j];
            hash = (hash << 2) | c;
        }
        kv_push(u64, *hash_list, hash);
    }
    return (int)kv_size(*hash_list);
}

void
collect_seeds(vec_u64* hash_list,
        const int read_id,
        const int read_start_id,
        const int reference_start_id,
        PackedDB* reference,
        LookupTable* lktbl,
        const int block_size,
		const int scan_window,
        const BOOL pairwise,
        vec_chain_seed* seed_list)
{
    idx soff_max = IDX_MAX;
    if (pairwise) {
        int max_rid = reference_start_id + (int)PDB_NUM_SEQS(reference);
        if (read_id + read_start_id >= reference_start_id && read_id + read_start_id < max_rid) {
            soff_max = PDB_SEQ_OFFSET(reference, read_id);
        }
    }
    const int n_kmer = (int)kv_size(*hash_list);
    ChainSeed seed;
    kv_clear(*seed_list);
    for (int i = 0; i < n_kmer; ++i) {
        seed.qoff = i * scan_window;
        u64 n_match = 0;
        u64* offset_list = extract_kmer_list(lktbl, kv_A(*hash_list, i), &n_match);
        if (offset_list) {
            for (u64 j = 0; j < n_match; ++j) {
                if (offset_list[j] >= soff_max) continue;
				seed.soff = offset_list[j];
                kv_push(ChainSeed, *seed_list, seed);
            }
        }
    }
}

void
build_block_seeds_index(vec_chain_seed* seed_list,
						const size_t block_size,
						vec_bseeds_index* index_list)
{
	kv_clear(*index_list);
	const size_t n_seeds = kv_size(*seed_list);
	ChainSeed* csv = kv_data(*seed_list);
	ks_introsort(ChainSeed, n_seeds, csv);
	
	size_t i = 0;
	BlockSeedsIndex bsidx;
	while (i < n_seeds) {
		size_t block_id = csv[i].soff / block_size; 
		size_t block_end = (block_id + 1) * block_size;
		size_t j =  i + 1;
		while (j < n_seeds && csv[j].soff < block_end) ++j;
		bsidx.start = i;
		bsidx.cnt = j - i;
		bsidx.block_id = block_id;
		kv_push(BlockSeedsIndex, *index_list, bsidx);
		i = j;
	}
}

int
find_max_score_seed(ChainSeed* seed_list,
					const size_t n_seeds,
					const double ddfs_cutoff,
					const int block_score_cutoff,
					vec_int* scores,
					size_t* chosen_id)
{
	kv_resize(int, *scores, n_seeds);
	memset(kv_data(*scores), 0, sizeof(int) * n_seeds);
	for (size_t i = 0; i < n_seeds - 1; ++i) {
		for (size_t j = i + 1; j < n_seeds; ++j) {
			if (ddf_score(seed_list[i].qoff, seed_list[i].soff, seed_list[j].qoff, seed_list[j].soff) < ddfs_cutoff) {
				++kv_A(*scores, i);
				++kv_A(*scores, j);
			}
		}
	}
	
	size_t max_i = 0;
	int max_score = kv_A(*scores, 0);
	for (size_t i = 1; i < n_seeds; ++i) {
		int score = kv_A(*scores, i);
		if (score > max_score) {
			max_i = i;
			score = max_score;
		}
	}
	if (max_score >= block_score_cutoff) {
		*chosen_id = max_i;
		return max_score;
	} else {
		*chosen_id = IDX_MAX;
		return 0;
	}
}

void
find_block_candidate(vec_chain_seed* seed_list,
					 vec_bseeds_index* index_list,
					 const size_t index_id,
					 PackedDB* reference,
					 const double ddfs_cutoff,
					 const int block_score_cutoff,
					 const size_t block_size,
					 vec_int* scores,
					 vec_chain_seed* chain_seed_list,
					 ChainDpData* chain_data,
					 const int qid,
					 const int qdir,
					 const idx qsize,
					 vec_can* candidates)
{
	size_t block_id = kv_A(*index_list, index_id).block_id;
	ChainSeed* block_seed_list = kv_data(*seed_list) + kv_A(*index_list, index_id).start;
	size_t max_score_id;
	if (!find_max_score_seed(block_seed_list, kv_A(*index_list, index_id).cnt, ddfs_cutoff, block_score_cutoff, scores, &max_score_id)) return;
	idx can_qoff = block_seed_list[max_score_id].qoff;
	idx can_soff = block_seed_list[max_score_id].soff;
	idx sid = pdb_offset_to_id(reference, can_soff);
	idx ssize = PDB_SEQ_SIZE(reference, sid);
	idx sref = PDB_SEQ_OFFSET(reference, sid);
	idx eref = sref + PDB_SEQ_SIZE(reference, sid);
	idx n = (can_soff - sref + block_size - 1) / block_size;
	idx sbid = (block_id >= n) ? (block_id - n) : 0;
	n = (eref - can_soff + block_size - 1) / block_size;
	idx ebid = block_id + n;
	
	idx index_id_from = index_id;
	while (1) {
		block_id = kv_A(*index_list, index_id_from).block_id;
		if (block_id < sbid) {
			++index_id_from;
			break;
		}
		if (index_id_from == 0) break;
		--index_id_from;
	}
	idx index_id_to = index_id;
	while (1) {
		block_id = kv_A(*index_list, index_id_to).block_id;
		if (block_id > ebid) break;
		if (index_id_to >= kv_size(*index_list)) break;
		++index_id_to;
	}
	
	kv_clear(*chain_seed_list);
	for (size_t i = index_id_from; i < index_id_to; ++i) {
		block_id = kv_A(*index_list, i).block_id;
		idx n_seeds = kv_A(*index_list, i).cnt;
		if (n_seeds == 0) continue;
		idx cnt = 0;
		block_seed_list = kv_data(*seed_list) + kv_A(*index_list, i).start;
		for (size_t k = 0; k < n_seeds; ++k) {
			if (block_seed_list[k].qoff == can_qoff && block_seed_list[k].soff == can_soff) {
				kv_push(ChainSeed, *chain_seed_list, block_seed_list[k]);
				++cnt;
				continue;
			}
			if (block_seed_list[k].soff < sref || block_seed_list[k].soff >= eref) continue;
			if (ddf_score(can_qoff, can_soff, block_seed_list[k].qoff, block_seed_list[k].soff) < ddfs_cutoff) {
				kv_push(ChainSeed, *chain_seed_list, block_seed_list[k]);
				++cnt;
			}
		}
		if (cnt * 10 >= n_seeds / 4) kv_A(*index_list, i).cnt = 0;
	}
	
	size_t n_chain_seeds = kv_size(*chain_seed_list);
	for (size_t i = 0; i < n_chain_seeds; ++i) {
		kv_A(*chain_seed_list, i).soff -= sref;
	}
	//OC_LOG("number of seeds: %lu", n_chain_seeds);
	ks_introsort(ChainSeed, n_chain_seeds, kv_data(*chain_seed_list));
	chain_data->chain_seeds = chain_seed_list;
	//for (size_t i = 0; i < n_chain_seeds; ++i) {
	//	OC_LOG("qoff = %lu, soff = %lu", kv_A(*chain_seed_list, i).qoff, kv_A(*chain_seed_list, i).soff);
	//}
	chain_dp(chain_data, qid, qdir, qsize, sid, ssize);
	//OC_LOG("number of candidates: %lu", kv_size(chain_data->lcanv));
	for (size_t i = 0; i != kv_size(chain_data->lcanv); ++i) {
		kv_push(GappedCandidate, *candidates, kv_A(chain_data->lcanv, i));
	}
}

WordFindData*
new_WordFindData(int kmer_size, int block_score_cutoff)
{
	WordFindData* data = (WordFindData*)malloc( sizeof(WordFindData) );
	kv_init(data->hash_list);
	kv_init(data->seed_list);
	kv_init(data->chain_seed_list);
	kv_init(data->ddf_scores);
	kv_init(data->candidates);
	kv_init(data->block_seeds_index);
	data->chain_data = new_ChainDpData(kmer_size, block_score_cutoff);
	
	return data;
}

WordFindData*
free_WordFindData(WordFindData* data)
{
	kv_destroy(data->hash_list);
	kv_destroy(data->seed_list);
	kv_destroy(data->chain_seed_list);
	kv_destroy(data->ddf_scores);
	kv_destroy(data->candidates);
	kv_destroy(data->block_seeds_index);
	data->chain_data = free_ChainDpData(data->chain_data);
	free(data);
	return 0;
}

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
			    vec_can* candidates)
{
	extract_hash_values(read, 
						read_size, 
						options->kmer_size, 
						options->scan_window, 
						&wfdata->hash_list);
	
	collect_seeds(&wfdata->hash_list, 
				  qid, 
				  read_start_id, 
				  reference_start_id, 
				  reference, 
				  lktbl, 
				  options->block_size, 
				  options->scan_window, 
				  pairwise, 
				  &wfdata->seed_list);
	
	build_block_seeds_index(&wfdata->seed_list, options->block_size, &wfdata->block_seeds_index);
	
	size_t n_blocks = kv_size(wfdata->block_seeds_index);
	//OC_LOG("number of blocks: %lu", n_blocks);
	for (size_t i = 0; i != n_blocks; ++i) {
		if (kv_A(wfdata->block_seeds_index, i).cnt > options->block_score_cutoff) {
			find_block_candidate(&wfdata->seed_list, 
								 &wfdata->block_seeds_index,
								 i,
								 reference,
								 options->ddfs_cutoff,
								 options->block_score_cutoff,
								 options->block_size,
								 &wfdata->ddf_scores,
								 &wfdata->chain_seed_list,
								 wfdata->chain_data,
								 qid,
								 qdir,
								 read_size,
								 candidates);
		}
	}
}
