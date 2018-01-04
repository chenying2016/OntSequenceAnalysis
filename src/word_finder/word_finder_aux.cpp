#include "word_finder_aux.h"

#include <algorithm>
#include <set>
#include <sstream>

using namespace std;

void 
GappedCandidates::filter_candidates(const int m) {
	int n = raw_cans.size();
	if (n <= m) {
		cans.assign(raw_cans.begin(), raw_cans.end());
		return;
	}
	
	sort(raw_cans.begin(), raw_cans.end(), GappedCandidateScore_GT());
	set<int> used_ids;
	int c = 0;
	for (size_t i = 0; i != raw_cans.size(); ++i) {
		int qid = raw_cans[i].qid;
		if (c < m) {
			cans.push_back(raw_cans[i]);
			if (used_ids.find(qid) == used_ids.end()) {
				++c;
				used_ids.insert(qid);
			}
		} else {
			if (used_ids.find(qid) != used_ids.end()) {
				cans.push_back(raw_cans[i]);
			}
		}
	}
}
