#include "results_pool.h"

#include <algorithm>

using namespace std;

const int EditScriptWorker::kMaxNum;

void
add_one_edit_segment(const u8 type, int num, std::vector<u8>& scripts)
{
	const int MaxSize = EditScriptWorker::kMaxNum;
	int left = num;
	while (left) {
		int n = min(left, MaxSize);
		u8 es = EditScriptWorker::pack_edit_script(type, n);
		scripts.push_back(es);
		left -= n;
	}
}

void 
build_edit_scripts(const std::string& qaln, const std::string& taln, std::vector<u8>& scripts)
{
	scripts.clear();
	r_assert(qaln.size() == taln.size())(qaln.size())(taln.size());
	int i = 0, j = 0, n = qaln.size();
	while (i < n) {
		if (qaln[i] == '-') {
			j = i + 1;
			while (j < n && qaln[j] == '-') ++j;
			add_one_edit_segment(EditScriptWorker::kOpD, j - i, scripts);
		} else if (taln[i] == '-') {
			j = i + 1;
			while (j < n && taln[j] == '-') ++j;
			add_one_edit_segment(EditScriptWorker::kOpI, j - i, scripts);
		} else {
			j = i + 1;
			while (j < n && qaln[j] != '-' && taln[j] != '-') ++j;
			add_one_edit_segment(EditScriptWorker::kOpM, j - i, scripts);
		}
		i = j;
	}
}

inline bool
interval_contained(idx ll, idx lr, idx sl, idx sr, int e)
{
	return (sl + e >= ll) && (sr <= lr + e);
}

void
ResultsPool::filter_contained_results()
{
	sort(results.begin(), results.end(), CmpMappingResult_Range(m4v));
	int n = results.size();
	for (int i = 0; i < n - 1; ++i) {
		if (!results[i].valid) continue;
		M4Record& mi = get_m4(results[i].m4_idx);
		for (int j = i + 1; j < n; ++j) {
			if (!results[j].valid) continue;
			M4Record& mj = get_m4(results[j].m4_idx);
			if (mi.sid != mj.sid) continue;
			if (mi.qdir != mj.qdir) continue;
			if (!interval_contained(mi.qoff, mi.qend, mj.qoff, mj.qend, 100)) continue;
			if (!interval_contained(mi.soff, mi.send, mj.soff, mj.send, 100)) continue;
			results[j].valid = 0;
		}
	}
	
	int i = 0, j = 0;
	for (; i < n; ++i) {
		if (results[i].valid) {
			if (i > j) results[j] = results[i];
			++j;
		}
	}
	
	results.resize(j);
}
