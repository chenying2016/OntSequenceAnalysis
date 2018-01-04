#ifndef RESULTS_POOL_H
#define RESULTS_POOL_H

#include "../common/m4record.h"
#include "../common/aux_tools.h"
#include "../common/pod_darr.h"
#include "../word_finder/word_finder_aux.h"

#include <vector>

struct EditScriptWorker
{
	static u8 pack_edit_script(const u8 type, const u8 num) {
		r_assert(num <= kMaxNum);
		return type | num;
	}
	static void unpack_edit_script(const u8 es, int& type, int& num) {
		type = es & kTypeMask;
		num = es & kNumMask;
	}
	
	static const int kMaxNum = 63;
	static const u8 kOpM = 64;
	static const u8 kOpI = 128;
	static const u8 kOpD = 192;
	static const u8 kTypeMask = 192;
	static const u8 kNumMask = 63;	
};

void
add_one_edit_segment(const u8 type, int num, std::vector<u8>& scripts);

void 
build_edit_scripts(const std::string& qaln, const std::string& taln, std::vector<u8>& scripts);

struct MappingResult
{
	int es_idx, es_size;
	int valid;
	int m4_idx, prev, next, parent;
};

#define can_contained(l, r, a, e) (((a) + (e) >= (l)) && ((a) <= (r) + (e)))

struct ResultsPool
{
	PODArray<MappingResult> results;
	PODArray<M4Record> m4v;
	PODArray<u8> edit_scripts;
	std::vector<u8> esv;
	
	ResultsPool() : edit_scripts(1000000) {}
	
	void clear() {
		results.clear();
		m4v.clear();
		edit_scripts.clear();
	}
	
	M4Record& get_m4(const int i) {
		return m4v[i];
	}
	
	bool check_candidate_contained(GappedCandidate& c) {
		for (int i = 0; i < results.size(); ++i) {
			M4Record& m = get_m4(results[i].m4_idx);
			if (m.sid != c.sid) continue;
			if (m.qdir != c. qdir) continue;
			if (!can_contained(m.qoff, m.qend, c.qoff, 100)) continue;
			if (!can_contained(m.soff, m.send, c.soff, 100)) continue;
			return true;
		}
		return false;
	}
	
	void filter_contained_results();
	
	void add_one_result(M4Record& m,
						const std::string& qaln,
						const std::string& saln) {
		MappingResult r;
		r.m4_idx = m4v.size();
		build_edit_scripts(qaln, saln, esv);
		r.es_idx = edit_scripts.size();
		r.es_size = esv.size();
		r.valid = 1;
		r.prev = -1;
		r.next = -1;
		r.parent = -1;
		m4v.push_back(m);
		results.push_back(r);
		edit_scripts.push_back(esv.data(), esv.size());
	}
	
	template <class W>
	void output(W& out, const int num_output) {
		if (results.empty()) return;
		int noutput = 0;
		int nrs = results.size();
		for (int i = 0; i < nrs && noutput < num_output; ++i) {
			MappingResult& r = results[i];
			if (r.parent != -1) continue;
			M4Record& mh = m4v[ r.m4_idx ];
			out << mh;
			if (r.prev != -1) {
				M4Record& mp = m4v[ r.prev ];
				out << mp;
			}
			if (r.next != -1) {
				M4Record& mn = m4v[ r.next ];
				out << mn;
			}
			++noutput;
		}
	}
};

struct CmpMappingResult_Range
{
	bool operator()(const MappingResult& a, const MappingResult& b) const {
		M4Record& m1 = m4v[a.m4_idx];
		M4Record& m2 = m4v[b.m4_idx];
		int r1 = m1.qend - m1.qoff;
		int r2 = m2.qend - m2.qoff;
		return r1 > r2;
	}
	
	CmpMappingResult_Range(PODArray<M4Record>& v) : m4v(v) {}
	
	PODArray<M4Record>& m4v;
};

#endif // RESULTS_POOL_H
