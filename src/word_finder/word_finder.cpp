#include "word_finder.h"

#include "chain_dp.h"

#include <algorithm>

using namespace std;

int
WordFinder::extract_hash_values()
{
	kmer_hash_values.clear();
    const int kmer_size = param->kmer_size;
    const int scan_stride = param->scan_stride;
    int from = 0, to = read_size - kmer_size;
    for (int i = from; i <= to; i += scan_stride) {
        u32 v = 0;
        for (int j = 0; j < kmer_size; ++j) {
            const u8 c = (u8)read[i + j];
            v = (v << 2) | c;
        }
        kmer_hash_values.push_back(v);
    }
    return kmer_hash_values.size();
}

void
print_kmer(u32 u, int kmer_size)
{
	string s;
	cout << u << '\t';
	for (int i = 0; i < kmer_size; ++i) {
		u8 c = u & 3;
		s += "ACGT"[c];
		u >>= 2;
	}
	for (size_t z = s.size(); z > 0; --z) {
		cout << s[z - 1];
	}
	cout << '\n';
}

void
WordFinder::collect_seeds(const bool pairwise)
{
    const int block_size = param->block_size;
    const int nkmer = extract_hash_values();
	idx soff_max = IDX_MAX;
	if (pairwise && read_id >= reference_start_id && read_id < reference_start_id + reference->num_seqs()) {
		soff_max = reference->seq_offset(read_id - reference_start_id);
	}
    SeedMatch sm;
    LookupTable::index_t* arr;
    int nsm;
    seeds.clear();
    for (int i = 0; i < nkmer; ++i) {
        sm.qno = i;
        if (lktbl->extract_kmer_list(kmer_hash_values[i], arr, nsm)) {
            for (int j = 0; j < nsm; ++j) {
				if (arr[j] >= soff_max) continue;
                sm.bid = static_cast<int>(arr[j] / block_size);
                sm.boff = arr[j] % block_size;
                seeds.push_back(sm);
            }
        }
    }
    sort(seeds.begin(), seeds.end(), CmpSeedMatch_Offset());
}

void 
WordFinder::calc_ddf_scores(const SeedMatch* smv, const int n)
{
    ddf_scores.resize(n);
    for (int i = 0; i < n; ++i) {
        ddf_scores[i].id = i;
        ddf_scores[i].score = 0;
    }
    int max_s_distance = static_cast<int>(read_size * (1.0 + param->error_rate));
    int scan_stride = param->scan_stride;
	double ddfs_cutoff = param->ddfs_cutoff;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
			if (check_ddf_score(smv[i].qno, smv[j].qno, smv[i].boff, smv[j].boff, scan_stride, max_s_distance, ddfs_cutoff)) {
				++ddf_scores[i].score;
				++ddf_scores[j].score;
			}
        }
    }
}

void
WordFinder::fill_one_block_seeds(const SeedMatch* smv, const int n)
{
    const int sm = SM;
    const int bid = smv[0].bid;
    BlockSeeds* block = read_seeds->block_seeds + bid;

    if (n <= sm) {
        for (int i = 0; i < n; ++i) {
            block->bseeds[i].qno = smv[i].qno;
            block->bseeds[i].boff = smv[i].boff;
        }
    } else {
        calc_ddf_scores(smv, n);
        sort(ddf_scores.begin(), ddf_scores.end(), CmpSeedDDFScore_Score());
        sort(ddf_scores.begin(), ddf_scores.begin() + sm, CmpSeedDDFScore_Id());
        for (int i = 0; i < sm; ++i) {
            int id = ddf_scores[i].id;
            block->bseeds[i].qno = smv[id].qno;
            block->bseeds[i].boff = smv[id].boff;
        }
    }

    read_seeds->block_ids[read_seeds->used_blocks] = bid;
    block->score = (short)n;
    block->score2 = (short)n;
    ++read_seeds->used_blocks;
}

void
WordFinder::fill_block_seeds()
{
    idx n = seeds.size();
    int i = 0;
    while (i < n) {
        int j = i + 1;
        while (j < n && seeds[i].bid == seeds[j].bid) ++j;
        fill_one_block_seeds(seeds.data() + i, j - i);
        i = j;
    }
}

int
find_gapped_start_locations(int* sids,
							int* boffs,
							int* scores,
							const int n,
							int scan_stride,
							int max_s_distance,
							int score_cutoff,
							double ddfs_cutoff)
{
	fill(scores, scores + n, 0);
	for (int i = 0; i < n - 1; ++i) 
		for (int j = i + 1, tempi = sids[i]; j < n; ++j) {
			if (tempi != sids[j] && check_ddf_score(sids[i], sids[j], boffs[i], boffs[j], scan_stride, max_s_distance, ddfs_cutoff)) {
				++scores[i];
				++scores[j];
				tempi = sids[j];
			}
		}
	
	int max_score = 0;
	int maxi = -1;
	int rep = 0;
	int lasti = -1;
	for (int i = 0; i < n; ++i)
		if (max_score < scores[i]) {
			max_score = scores[i];
			maxi = i;
			rep = 0;
		} else if (max_score == scores[i]) {
			++rep;
			lasti = i;
		}
	if (max_score < score_cutoff) return -1;
	
	if (rep == max_score) return maxi;
	
	for (int i = 0; i < maxi; ++i)
		if (check_ddf_score(sids[i], sids[maxi], boffs[i], boffs[maxi], scan_stride, max_s_distance, ddfs_cutoff)) {
			return i;
		}
	
	for (int i = maxi + 1; i < n; ++i)
		if (check_ddf_score(sids[maxi], sids[i], boffs[maxi], boffs[i], scan_stride, max_s_distance, ddfs_cutoff)) {
			return i;
		}
	
	return -1;
}

inline bool
check_can_range_ratio(GappedCandidate& can, double ratio)
{
	idx L1 = can.qoff, R1 = can.qsize - can.qoff;
	idx L2 = can.soff, R2 = can.ssize - can.soff;
	idx L = min(L1, L2);
	idx R = min(R1, R2);
	idx S = L + R;
	return (S >= can.qsize * ratio) || (S >= can.ssize * ratio);
}

void 
WordFinder::calc_candidate_for_one_block(const int block_id, GappedCandidates& cans)
{
	const int sm = SM;
	const int block_size = param->block_size;
	const int scan_stride = param->scan_stride;
	const int score_cutoff = param->block_score_cutoff;
	const double ddfs_cutoff = param->ddfs_cutoff;
	int sids[sm * 2], boffs[sm * 2], scores[sm * 2];
	BlockSeeds* block = read_seeds->block_seeds + block_id;
	int nseeds = 0;
	for (int i = 0; i < block->score && i < sm; ++i) {
		sids[nseeds] = block->bseeds[i].qno;
		r_assert(sids[nseeds] * scan_stride < read_size);
		boffs[nseeds] = block->bseeds[i].boff;
		++nseeds;
	}
	++block;
	if (block->score) {
		for (int i = 0; i < block->score && i < sm; ++i) {
			sids[nseeds] = block->bseeds[i].qno;
			r_assert(sids[nseeds] * scan_stride < read_size);
			boffs[nseeds] = block->bseeds[i].boff + block_size;
			++nseeds;
		}
	}
	--block;
	
	const int max_s_distance = read_size * (1.0 + param->error_rate);
	int mid = find_gapped_start_locations(sids, boffs, scores, nseeds, scan_stride, max_s_distance, score_cutoff, ddfs_cutoff);
	if (mid == -1) return;
	
	idx tsoff = block_id; tsoff *= block_size; tsoff += boffs[mid];
	idx sid = reference->offset_to_id(tsoff);
	idx sstart = reference->seq_offset(sid);
	idx ssize = reference->seq_size(sid);
	idx send = sstart + ssize;
	
	int tsid = sids[mid];
	idx tqoff = sids[mid] * scan_stride;
	r_assert(tqoff < read_size);
	idx num1 = min(tqoff, tsoff);
	idx num2 = min(read_size - tqoff, reference->size() - tsoff);
	int score_add = 0;
	
	block = read_seeds->block_seeds + block_id - 1;
	idx n = (num1 + block_size - 1) / block_size;
	idx bid = block_id - 1;
	for (; n >= 0 && bid >= 0; --n, --bid, --block) {
		if (block->score) {
			idx block_start = bid * block_size;
			int scnt = min((int)block->score, sm);
			int qn = 0;
			for (int i = 0; i < scnt; ++i) {
				idx sdist = tsoff - block_start - block->bseeds[i].boff;
				double s = 1.0;
				if (sdist < max_s_distance && tsid != block->bseeds[i].qno) {
					s = 1.0 * sdist / ((tsid - block->bseeds[i].qno) * scan_stride);
					s = fabs(1.0 - s);
				}
				if (s < ddfs_cutoff) {
					++qn;
					++score_add;
				}
			}
			if (1.0 * qn / scnt >= 0.4) block->score = 0;
		}
	}
	
	block = read_seeds->block_seeds + block_id + 1;
	block->score = 0;
	++block;
	n = (num2 + block_size - 1) / block_size;
	bid = block_id + 2;
	for (; n >= 0; --n, ++bid, ++block) {
		if (block->score) {
			idx block_start = bid * block_size;
			int scnt = min((int)block->score, sm);
			int qn = 0;
			for (int i = 0; i < scnt; ++i) {
				idx sdist = block_start + block->bseeds[i].boff - tsoff;
				double s = 1.0;
				if (sdist < max_s_distance && tsid != block->bseeds[i].qno) {
					s = 1.0 * sdist / ((block->bseeds[i].qno - tsid) * scan_stride);
					s = fabs(1.0 - s);
				}
				if (s < ddfs_cutoff) {
					++qn;
					++score_add;
				}
			}
			if (1.0 * qn / scnt >= 0.4) block->score = 0;
		}
	}
	
	GappedCandidate can;
	can.qdir = read_dir;
	can.qoff = tqoff;
	can.qsize = read_size;
	can.score = scores[mid] + score_add;
	can.sid = sid;
	can.soff = tsoff - sstart;
	can.ssize = ssize;
	if (can.score < param->block_score_cutoff) return;
	
	if (1)
	{
		vector<ChainSeed>& cseeds = chain_data.seeds;
		cseeds.clear();
		n = (num1 + block_size - 1) / block_size;
		idx sbid = block_id - n;
		if (sbid < 0) sbid = 0;
		n = (num2 + block_size - 1) / block_size;
		idx ebid = block_id + n;
		idx sref = reference->seq_offset(sid);
		idx eref = reference->seq_size(sid) + sref;
		for (idx i = sbid; i <= ebid; ++i) {
			idx sblk = i * block_size;
			block = read_seeds->block_seeds + i;
			int scnt = min((int)block->score2, sm);
			for (int k = 0; k < scnt; ++k) {
				idx soff = sblk + block->bseeds[k].boff;
				if (soff < sref || soff >= eref) continue;
				soff -= sref;
				idx qoff = block->bseeds[k].qno;
				qoff *= scan_stride;
				if (qoff == can.qoff && soff == can.soff) {
					cseeds.push_back(ChainSeed(qoff, soff));
				} else if (soff == can.soff || qoff == can.qoff) {
					continue;
				}
				if ((qoff < can.qoff && soff < can.soff) || (qoff > can.qoff && soff > can.soff)) {
					double s = 1.0 * (soff - can.soff) / (qoff - can.qoff);
					s = fabs(s);
					s = fabs(1.0 - s);
					if (s <= ddfs_cutoff) {
						cseeds.push_back(ChainSeed(qoff, soff));
					}
				}
			}
		}
		
		sort(cseeds.begin(), cseeds.end(), ChainSeedSoffLt());
		chain_dp(chain_data, read_id, read_dir, read_size, sid, ssize);
		for (size_t z = 0; z != chain_data.lcanv.size(); ++z) {
			cans.add(chain_data.lcanv[z]);
		}
	}
}

void 
WordFinder::calc_candidate_for_one_block(const int block_id, CandidateStore& cans)
{
	const int sm = SM;
	const int block_size = param->block_size;
	const int scan_stride = param->scan_stride;
	const int score_cutoff = param->block_score_cutoff;
	const double ddfs_cutoff = param->ddfs_cutoff;
	int sids[sm * 2], boffs[sm * 2], scores[sm * 2];
	BlockSeeds* block = read_seeds->block_seeds + block_id;
	int nseeds = 0;
	for (int i = 0; i < block->score && i < sm; ++i) {
		sids[nseeds] = block->bseeds[i].qno;
		r_assert(sids[nseeds] * scan_stride < read_size);
		boffs[nseeds] = block->bseeds[i].boff;
		++nseeds;
	}
	++block;
	if (block->score) {
		for (int i = 0; i < block->score && i < sm; ++i) {
			sids[nseeds] = block->bseeds[i].qno;
			r_assert(sids[nseeds] * scan_stride < read_size);
			boffs[nseeds] = block->bseeds[i].boff + block_size;
			++nseeds;
		}
	}
	--block;
	
	const int max_s_distance = read_size * (1.0 + param->error_rate);
	int mid = find_gapped_start_locations(sids, boffs, scores, nseeds, scan_stride, max_s_distance, score_cutoff, ddfs_cutoff);
	if (mid == -1) return;
	
	idx tsoff = block_id; tsoff *= block_size; tsoff += boffs[mid];
	idx sid = reference->offset_to_id(tsoff);
	idx sstart = reference->seq_offset(sid);
	idx ssize = reference->seq_size(sid);
	idx send = sstart + ssize;
	
	int tsid = sids[mid];
	idx tqoff = sids[mid] * scan_stride;
	r_assert(tqoff < read_size);
	idx num1 = min(tqoff, tsoff);
	idx num2 = min(read_size - tqoff, reference->size() - tsoff);
	int score_add = 0;
	
	block = read_seeds->block_seeds + block_id - 1;
	idx n = (num1 + block_size - 1) / block_size;
	idx bid = block_id - 1;
	for (; n >= 0 && bid >= 0; --n, --bid, --block) {
		if (block->score) {
			idx block_start = bid * block_size;
			int scnt = min((int)block->score, sm);
			int qn = 0;
			for (int i = 0; i < scnt; ++i) {
				idx sdist = tsoff - block_start - block->bseeds[i].boff;
				double s = 1.0;
				if (sdist < max_s_distance && tsid != block->bseeds[i].qno) {
					s = 1.0 * sdist / ((tsid - block->bseeds[i].qno) * scan_stride);
					s = fabs(1.0 - s);
				}
				if (s < ddfs_cutoff) {
					++qn;
					++score_add;
				}
			}
			if (1.0 * qn / scnt >= 0.4) block->score = 0;
		}
	}
	
	block = read_seeds->block_seeds + block_id + 1;
	block->score = 0;
	++block;
	n = (num2 + block_size - 1) / block_size;
	bid = block_id + 2;
	for (; n >= 0; --n, ++bid, ++block) {
		if (block->score) {
			idx block_start = bid * block_size;
			int scnt = min((int)block->score, sm);
			int qn = 0;
			for (int i = 0; i < scnt; ++i) {
				idx sdist = block_start + block->bseeds[i].boff - tsoff;
				double s = 1.0;
				if (sdist < max_s_distance && tsid != block->bseeds[i].qno) {
					s = 1.0 * sdist / ((block->bseeds[i].qno - tsid) * scan_stride);
					s = fabs(1.0 - s);
				}
				if (s < ddfs_cutoff) {
					++qn;
					++score_add;
				}
			}
			if (1.0 * qn / scnt >= 0.4) block->score = 0;
		}
	}
	
	GappedCandidate can;
	can.qdir = read_dir;
	can.qoff = tqoff;
	can.qsize = read_size;
	can.score = scores[mid] + score_add;
	can.sid = sid;
	can.soff = tsoff - sstart;
	can.ssize = ssize;
	if (can.score < param->block_score_cutoff) return;
	
	if (1)
	{
		vector<ChainSeed>& cseeds = chain_data.seeds;
		cseeds.clear();
		n = (num1 + block_size - 1) / block_size;
		idx sbid = block_id - n;
		if (sbid < 0) sbid = 0;
		n = (num2 + block_size - 1) / block_size;
		idx ebid = block_id + n;
		idx sref = reference->seq_offset(sid);
		idx eref = reference->seq_size(sid) + sref;
		for (idx i = sbid; i <= ebid; ++i) {
			idx sblk = i * block_size;
			block = read_seeds->block_seeds + i;
			int scnt = min((int)block->score2, sm);
			for (int k = 0; k < scnt; ++k) {
				idx soff = sblk + block->bseeds[k].boff;
				if (soff < sref || soff >= eref) continue;
				soff -= sref;
				idx qoff = block->bseeds[k].qno;
				qoff *= scan_stride;
				if (qoff == can.qoff && soff == can.soff) {
					cseeds.push_back(ChainSeed(qoff, soff));
				} else if (soff == can.soff || qoff == can.qoff) {
					continue;
				}
				if ((qoff < can.qoff && soff < can.soff) || (qoff > can.qoff && soff > can.soff)) {
					double s = 1.0 * (soff - can.soff) / (qoff - can.qoff);
					s = fabs(s);
					s = fabs(1.0 - s);
					if (s <= ddfs_cutoff) {
						cseeds.push_back(ChainSeed(qoff, soff));
					}
				}
			}
		}
		
		sort(cseeds.begin(), cseeds.end(), ChainSeedSoffLt());
		chain_dp(chain_data, read_id, read_dir, read_size, sid, ssize);
		for (size_t z = 0; z != chain_data.lcanv.size(); ++z) {
			cans.add(chain_data.lcanv[z]);
		}
	}
}

void 
WordFinder::calc_candidates(GappedCandidates& cans)
{
    const int bscore_cutoff = param->block_score_cutoff;
    for (int i = 0; i < read_seeds->used_blocks; ++i) {
        int block_id = read_seeds->block_ids[i];
        if (read_seeds->block_seeds[block_id].score >= bscore_cutoff) {
            calc_candidate_for_one_block(block_id, cans);
        }
    }
}

void 
WordFinder::calc_candidates(CandidateStore& cans)
{
    const int bscore_cutoff = param->block_score_cutoff;
    for (int i = 0; i < read_seeds->used_blocks; ++i) {
        int block_id = read_seeds->block_ids[i];
        if (read_seeds->block_seeds[block_id].score >= bscore_cutoff) {
            calc_candidate_for_one_block(block_id, cans);
        }
    }
}
