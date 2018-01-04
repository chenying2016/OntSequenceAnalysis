#ifndef CHAIN_DP_H
#define CHAIN_DP_H

#include "../common/gapped_candidate.h"
#include "word_finder_aux.h"
#include <vector>

using std::vector;

typedef std::pair<int, int> PInt;

struct ChainDpData
{
	vector<int> f;
	vector<int> p;
	vector<int> t;
	vector<int> v;
	vector<PInt> u;
	vector<GappedCandidate> lcanv;
	vector<ChainSeed> seeds;
	int kmer_size;
	int max_dist;
	int bw;
	int max_skip;
	int min_cnt;
	int min_sc;
};

void
chain_dp(ChainDpData& data, const int qid, const int qdir, const idx qsize, const idx sid, const idx ssize);
	

#endif // CHAIN_DP_H
