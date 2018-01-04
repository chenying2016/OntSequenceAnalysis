#ifndef WORD_FINDER_AUX_H
#define WORD_FINDER_AUX_H

#include "../common/gapped_candidate.h"
#include "../common/lookup_table.h"
#include "../common/packed_db.h"

#include <cmath>

#define SM 40
#define SI (SM + 1)

struct InitWordParameters
{
	int kmer_size;
	int scan_stride;
	int num_candidates;
	int block_score_cutoff;
	int block_size;
	double error_rate;
	double mapped_ratio_cutoff;
	double ddfs_cutoff;
};

struct SeedMatch
{
    int qno;
    int bid;
    int boff;
};

struct CmpSeedMatch_Offset
{
    bool operator()(const SeedMatch& a, const SeedMatch& b) const {
        if (a.bid != b.bid) return a.bid < b.bid;
        if (a.qno != b.qno) return a.qno < b.qno;
        return a.boff < b.boff;
    }
};

struct ChainSeed
{
	idx qoff;
	idx soff;
	
	ChainSeed(idx a, idx b): qoff(a), soff(b) {}
};

struct ChainSeedSoffLt
{
	bool operator()(const ChainSeed& a, const ChainSeed& b) {
		return (a.soff == b.soff) ? (a.qoff < b.qoff) : (a.soff < b.soff);
	}
};

struct BlockSeedMatch
{
    int qno;
    int boff;
};

struct BlockSeeds
{
    BlockSeedMatch bseeds[SM];
    short score;
    short score2;
};

struct Seeds
{
    BlockSeeds* block_seeds;
    int* block_ids;
    int used_blocks;

    Seeds(const idx reference_size, const idx block_size) {
        const idx num_blocks = (reference_size + block_size - 1) / block_size + 5;
        snew(block_ids, int, num_blocks);
        snew(block_seeds, BlockSeeds, num_blocks);
		for (idx i = 0; i < num_blocks; ++i) {
			block_seeds[i].score = 0;
			block_seeds[i].score2 = 0;
		}
        used_blocks = 0;
    }

    ~Seeds() {
        sfree(block_ids);
        sfree(block_seeds);
    }

    void clear() {
        for (int i = 0; i < used_blocks; ++i) {
            int bid = block_ids[i];
            block_seeds[bid].score = 0;
            block_seeds[bid].score2 = 0;
        }
        used_blocks = 0;
    }
};

struct SeedDDFScore
{
    int id;
    int score;
};

struct CmpSeedDDFScore_Score
{
    bool operator()(const SeedDDFScore& a, const SeedDDFScore& b) const {
        return a.score > b.score;
    }
};

struct CmpSeedDDFScore_Id
{
    bool operator()(const SeedDDFScore& a, const SeedDDFScore& b) const {
        return a.id < b.id;
    }
};

struct GappedCandidates
{
	/*
	void add(const GappedCandidate& can) {
        //std::cout << can;
        if (n < m) {
            cans[n++] = can;
        } else {
            r_assert(n == m);
            int min_i = 0;
            int min_score = cans[0].score;
            for (int i = 1; i < n; ++i) {
                if (min_score > cans[i].score) {
                    min_i = i;
                    min_score = cans[i].score;
                }
            }
            if (min_score < can.score) {
                cans[min_i] = can;
            }
        }
	}
	int n, m;
	GappedCandidate* cans;
	*/
	
	void clear() {
		raw_cans.clear();
		cans.clear();
	}
	
	void add(const GappedCandidate& can) {
		raw_cans.push_back(can);
	}
	
	void filter_candidates(const int m);
	
	std::vector<GappedCandidate> raw_cans;
	std::vector<GappedCandidate> cans;
};

inline double
calc_ddf_score(int qoffi, int qoffj, int soffi, int soffj)
{
	double s = 1.0 * (soffj - soffi) / (qoffj - qoffi);
	s = fabs(1.0 - s);
	return s;
}

inline bool
check_ddf_score(int qnoi, int qnoj, int soffi, int soffj, int scan_stride, int max_s_distance, double ddfs_cutoff)
{
	if (qnoi >= qnoj
		||
		soffi >= soffj
		||
		soffj - soffi >= max_s_distance) return false;
	
	int qoffi = qnoi * scan_stride;
	int qoffj = qnoj * scan_stride;
	double s = calc_ddf_score(qoffi, qoffj, soffi, soffj);
	return s < ddfs_cutoff;
}

#endif // WORD_FINDER_AUX_H
