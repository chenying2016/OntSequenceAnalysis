#ifndef WORD_FINDER_H
#define WORD_FINDER_H

#include "candidate_store.h"
#include "word_finder_aux.h"
#include "chain_dp.h"

int
find_gapped_start_locations(int* sids,
int* boffs,
int* scores,
const int n,
int scan_stride,
int max_s_distance,
int score_cutoff,
double ddfs_cutoff);

class WordFinder
{
public:
	
	WordFinder(InitWordParameters* _param,
			  LookupTable* _lktbl,
			  PackedDB*	_reference,
			  int _reference_start_id) {
		param = _param;
		lktbl = _lktbl;
		reference = _reference;
		reference_start_id = _reference_start_id;
		
		read = NULL;
		kmer_hash_values.reserve(100000);
		seeds.reserve(100000);
		ddf_scores.reserve(100000);
		the_read_seeds[0] = new Seeds(reference->size(), param->block_size);
		the_read_seeds[1] = new Seeds(reference->size(), param->block_size);
		read_seeds = NULL;
		
		chain_data.kmer_size = param->kmer_size;
		chain_data.max_dist = 10000;
		chain_data.bw = 500;
		chain_data.max_skip = 25;
		chain_data.min_cnt = param->block_score_cutoff;
		chain_data.min_sc = 100;
	}
	
	~WordFinder() {
		delete the_read_seeds[0];
		delete the_read_seeds[1];
	}
	
	void clear() {
        seeds.clear();
        ddf_scores.clear();
        the_read_seeds[0]->clear();
        the_read_seeds[1]->clear();
        read_seeds = NULL;
    }
	
	void calc_candidates(const char* new_read, 
						 const int new_read_id, 
						 const int new_read_dir, 
						 const int new_read_size,
						 GappedCandidates& cans,
						 const bool pairwise) {
        clear();
        set_read(new_read, new_read_id, new_read_dir, new_read_size);
        extract_hash_values();
        collect_seeds(pairwise);
        fill_block_seeds();
        calc_candidates(cans);
    }
	
	void calc_candidates(const char* new_read, 
						 const int new_read_id, 
						 const int new_read_dir, 
						 const int new_read_size,
						 CandidateStore& cans,
						 const bool pairwise) {
        clear();
        set_read(new_read, new_read_id, new_read_dir, new_read_size);
        extract_hash_values();
        collect_seeds(pairwise);
        fill_block_seeds();
        calc_candidates(cans);
    }
	
	BlockSeeds* strand_seeds(int dir) {
		return the_read_seeds[dir]->block_seeds;
	}
	
private:
	/// find seeds
    int extract_hash_values();
    void collect_seeds(const bool pairwise);
    void calc_ddf_scores(const SeedMatch* smv, const int n);
    void fill_one_block_seeds(const SeedMatch* smv, const int n);
    void fill_block_seeds();
	
	/// find candidates
    void calc_candidate_for_one_block(const int block_id, GappedCandidates& cans);
	void calc_candidate_for_one_block(const int block_id, CandidateStore& cans);
    void calc_candidates(GappedCandidates& cans);
	void calc_candidates(CandidateStore& cans);
	
	void set_read(const char* new_read,
				  const int new_read_id,
				  const int new_read_dir,
				  const int new_read_size) {
        read = new_read;
        read_id = new_read_id;
        read_dir = new_read_dir;
        read_size = new_read_size;
        read_seeds = the_read_seeds[new_read_dir];
    }
	
	WordFinder(const WordFinder&);
	WordFinder& operator = (const WordFinder&);
	
private:
	InitWordParameters*			param;
	const char*					read;
	int							read_id;
	int							read_dir;
	int 						read_size;
	LookupTable*				lktbl;
	PackedDB*					reference;
	int							reference_start_id;
	PODArray<u32>				kmer_hash_values;
	PODArray<SeedMatch>			seeds;
	PODArray<SeedDDFScore>		ddf_scores;
	Seeds*						the_read_seeds[2];
	Seeds*						read_seeds;
	ChainDpData 				chain_data;
};

#endif // WORD_FINDER_H
