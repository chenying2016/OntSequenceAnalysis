#ifndef CANDIDATE_STORE_H
#define CANDIDATE_STORE_H

#include "../common/gapped_candidate.h"
#include "../common/smart_assert.h"

#include <algorithm>
#include <set>
#include <vector>

using std::set;
using std::vector;

struct CandidateId
{
    std::size_t id;
    int score;
};

struct CandidateId_LT
{
    bool operator()(const CandidateId& a, const CandidateId& b) const {
        return a.score > b.score;
    }
};

struct CandidateStore
{
public:
	CandidateStore(size_t m): max_num_good_cans(m), first_time_full(1) {}
	
    void clear() {
        first_time_full = 1;
        canv.clear();
		good_canv.clear();
        good_can_ids.clear();
    }

    void add(const GappedCandidate& can) {
        CandidateId ci; ci.id = canv.size(); ci.score = can.score;
        canv.push_back(can);

        if (good_can_ids.size() != max_num_good_cans) {
            good_can_ids.push_back(ci);
        } else {
            r_assert(good_can_ids.size() == max_num_good_cans);
            if (first_time_full) {
                first_time_full = 0;
                make_heap(good_can_ids.begin(), good_can_ids.end(), CandidateId_LT());                
                for (size_t i = 1; i < max_num_good_cans; ++i) {
                    r_assert(good_can_ids[i].score >= good_can_ids[0].score);
                }
            }
            if (good_can_ids[0].score < can.score) {
                pop_heap(good_can_ids.begin(), good_can_ids.end(), CandidateId_LT());
                r_assert(good_can_ids[max_num_good_cans - 1].score < can.score);
                good_can_ids[max_num_good_cans - 1] = ci;
                push_heap(good_can_ids.begin(), good_can_ids.end(), CandidateId_LT());
            }
        }
    }

    template <class OutputStream>
        void output(OutputStream& out) {
            for (size_t i = 0; i != good_can_ids.size(); ++i) {
                GappedCandidate& can = canv[ good_can_ids[i].id ];
                out << can;
                can.score = 0;
            }
            if (canv.size() > max_num_good_cans) {
                set<int> qids;
                for (size_t i = 0; i != good_can_ids.size(); ++i) {
                    GappedCandidate& can = canv[ good_can_ids[i].id ];
                    qids.insert(can.qid);
                }
                for (size_t i = 0; i != canv.size(); ++i) {
                    if (canv[i].score == 0) continue;
                    if (qids.find(canv[i].score) != qids.end()) {
                        out << canv[i];
                    }
                }
            }
        }

    void fix_sequence_ids(const idx qs_id, const idx ss_id) {
        for (size_t i = 0; i != canv.size(); ++i) {
            canv[i].qid += qs_id;
            canv[i].sid += ss_id;
        }
    }
	
	vector<GappedCandidate>&
	get_good_candidates() {
		for (size_t i = 0; i != good_can_ids.size(); ++i) {
			good_canv.push_back(canv[good_can_ids[i].id]);
		}
		return good_canv;
	}

private:
    vector<GappedCandidate> canv;
	vector<GappedCandidate> good_canv;
    vector<CandidateId> good_can_ids;
    const std::size_t max_num_good_cans;
    bool first_time_full;
};

#endif // CANDIDATE_STORE_H
