#include "chain_dp.h"

#include "word_finder_aux.h"
#include "word_finder.h"

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
	-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
	LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
	LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
	uint32_t t, tt;
	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}

struct CmpPInt
{
    bool operator()(const PInt& a, const PInt& b) {
        return (a.first == b.first) ? (a.second < b.second) : (a.first > b.first);
    }
};

void
chain_dp(ChainDpData& data, const int qid, const int qdir, const idx qsize, const idx sid, const idx ssize)
{
	vector<int>& f = data.f;
	vector<int>& p = data.p;
	vector<int>& t = data.t;
	vector<int>& v = data.v;
	vector<PInt>& u = data.u;
	vector<GappedCandidate>& lcanv = data.lcanv;
	vector<ChainSeed>& seeds = data.seeds;
	lcanv.clear();
	const int kmer_size = data.kmer_size;
	const int max_dist = data.max_dist;
	const int bw = data.bw;
	const int max_skip = data.max_skip;
	const int min_cnt = data.min_cnt;
	const int min_sc = data.min_sc;
    const int n = seeds.size();
	
    f.resize(n); fill(f.begin(), f.end(), 0);
    p.resize(n); fill(p.begin(), p.end(), -1);
    t.resize(n); fill(t.begin(), t.end(), 0);
    v.resize(n); fill(v.begin(), v.end(), 0);
    int i, j, k, st = 0;

    for (i = 0; i < n; ++i) {
        idx ri = seeds[i].soff;
        int max_j = -1;
        idx qi = seeds[i].qoff;
        int max_f = kmer_size;
        int n_skip = 0;
        while (st < i && ri - seeds[st].soff > max_dist) ++st;
        for (j = i - 1; j >= st; --j) {
            idx dr = ri - seeds[j].soff;
            idx dq = qi - seeds[j].qoff;
            if (dq > max_dist) continue;
            const idx dd = dr > dq ? dr - dq : dq - dr;
            if (dd > bw) continue;
            const idx min_d = min(dq, dr);
            int sc = min_d;
            int log_dd = dd ? ilog2_32(dd) : 0;
            sc -= (int)(dd * 0.01 * kmer_size) + (log_dd >> 1);
            sc += f[j];
            if (sc > max_f) {
                max_f = sc;
                max_j = j;
                if (n_skip > 0) --n_skip;
            } else if (t[j] == i) {
                if (++n_skip > max_skip) break;
            }
            if (p[j] >= 0) t[p[j]] = i;
        }
        f[i] = max_f;
        p[i] = max_j;
        v[i] = max_j >= 0 && v[max_j] > max_f ? v[max_j] : max_f;
    }

    fill(t.begin(), t.end(), 0);
    for (i = 0; i < n; ++i) 
        if (p[i] >= 0) t[p[i]] = 1;
    int n_u;
    for (i = n_u = 0; i < n; ++i) {
        if (t[i] == 0 && v[i] >= min_sc) ++n_u;
    }
    if (n_u == 0) return;
    
    u.clear();
    for (i = n_u = 0; i < n; ++i) {
        if (t[i] == 0 && v[i] >= min_sc) {
            j = i;
            while (j >= 0 && f[j] < v[j]) j = p[j];
            if (j < 0) j = i;
            u.push_back(make_pair(f[j], j));
            ++n_u;
        }
    }
    sort(u.begin(), u.end(), CmpPInt());

	GappedCandidate can;
	can.qid = qid;
	can.qdir = qdir;
	can.qsize = qsize;
	can.sid = sid;
	can.sdir = FWD;
	can.ssize = ssize;

    fill(t.begin(), t.end(), 0);
    int n_v;
    for (i = n_v = k = 0; i < n_u; ++i) {
        int n_v0 = n_v, k0 = k;
        j = u[i].second;
		//can.qoff = seeds[j].qoff;
		//can.soff = seeds[j].soff;
		can.qend = seeds[j].qoff + kmer_size;
		can.send = seeds[j].soff + kmer_size;
		can.qoff = can.qend;
		can.soff = can.send;
        do {
            n_v++;
            t[j] = 1;
            j = p[j];
        } while (j >= 0 && t[j] == 0);
        if (j < 0) {
            if (n_v - n_v0 >= min_cnt) {
				can.qbeg = seeds[0].qoff;
				can.sbeg = seeds[0].soff;
				can.score = u[i].first;
				lcanv.push_back(can);
			}
        } else if (u[i].first - f[j] >= min_sc) {
            if (n_v - n_v0 >= min_cnt) {
				can.qbeg = seeds[n_v0].qoff;
				can.sbeg = seeds[n_v0].soff;
				can.score = u[i].first - f[j];
				lcanv.push_back(can);
			}
        }
        if (k0 == k) n_v = n_v0;
    }
	if (lcanv.empty()) return;
	
	sort(lcanv.begin(), lcanv.end(), GappedCandidateScore_GT());
	const int m = lcanv.size();
	for (i = 0; i < m; ++i) {
		if (lcanv[i].qid == -1) continue;
		for (j = i + 1; j < m; ++j) {
			if (lcanv[j].qid == -1) continue;
			if (lcanv[j].qoff >= lcanv[i].qbeg && lcanv[j].qoff <= lcanv[i].qend) {
				lcanv[j].qid = -1;
				continue;
			}
			if (lcanv[j].soff >= lcanv[i].sbeg && lcanv[j].soff <= lcanv[i].send) {
				lcanv[j].qid = -1;
				continue;
			}
		}
	}
	i = 0;
	for (j = 0; j < m; ++j) {
		if (lcanv[j].qid != -1) lcanv[i++] = lcanv[j];
	}
	lcanv.resize(i);
}
