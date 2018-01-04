#include "cns_stage1.h"

#include "../common/cns_seq.h"
#include "../common/smart_assert.h"
#include "consensus_aux.h"

#include <algorithm>
#include <set>
#include <string>

using namespace std;

struct Stage1CovCnt
{
	int cnts[6];
	
	int cov() const {
		return cnts[5];
	}
	
	Stage1CovCnt() {
		std::fill(cnts, cnts + 6, 0);
	}
};

void
add_one_align(const string& qaln, const string& taln, int toff, vector<Stage1CovCnt>& cnts)
{
	r_assert(qaln.size() == taln.size());
	for (size_t i = 0; i != qaln.size(); ++i) {
		if (taln[i] == GAP_CHAR) continue;
		++cnts[toff].cnts[5];
		switch (qaln[i]) {
			case GAP_CHAR:
				++cnts[toff++].cnts[4];
				break;
			case 'A':
				++cnts[toff++].cnts[0];
				break;
			case 'C':
				++cnts[toff++].cnts[1];
				break;
			case 'G':
				++cnts[toff++].cnts[2];
				break;
			case 'T':
				++cnts[toff++].cnts[3];
				break;
			default:
				break;
		}
	}
}

void
process_one_frag(vector<Stage1CovCnt>& cnts,
				 int read_id,
				 int read_size,
				 int from,
				 int to,
				 vector<char>& cns_seq)
{
	for (size_t i = from; i < to; ++i) {
		int max_cnt = 0, max_k = -1;
		for (int k = 0; k < 5; ++k) {
			if (cnts[i].cnts[k] > max_cnt) {
				max_cnt = cnts[i].cnts[k];
				max_k = k;
			}
		}
		if (max_k >= 0 && max_cnt >= cnts[i].cnts[5] * 0.6) {
			cns_seq[i] = max_k;
		}
	}
}

void
cns_stage1(ConsensusData* cns_data,
		   const int template_id,
		   const int cani_sid,
		   const int cani_eid,
		   const double min_mapping_ratio)
{
	OutputStream& out = cns_data->out;
    EGCInfo* caniv = cns_data->caniv;
    const int min_align_size = cns_data->options.min_align_size;
	const int min_cov = cns_data->options.min_cov;
	const int min_size = cns_data->options.min_size;
	OverlapsProvider& op = cns_data->op;
	set<int>& corrected_read_ids = cns_data->corrected_read_ids;
	const int template_size = caniv[cani_sid].template_size();
    int num_added = 0;
    int num_ext = 0;
    const int max_ext = MAX_CNS_OVLPS + 100;
    set<int> used_ids;
	int qid, qoff, qend, toff, tend;

	op.prefetch_target_sequences(template_id);
    sort(caniv + cani_sid, caniv + cani_eid, CmpEGCInfo_Score());
	
	vector<Stage1CovCnt> cov_cnts(template_size);
    for (int i = cani_sid; i < cani_eid && num_added < MAX_CNS_OVLPS && num_ext < max_ext; ++i, ++num_ext) {
		qid = caniv[i].non_template_id();
		if (used_ids.find(qid) != used_ids.end()) continue;
        if (op.get_overlaps(caniv[i], min_align_size, min_mapping_ratio, qid, qoff, qend, toff, tend, 0)) {
			if (!check_overlap_relations(qoff, qend, caniv[i].non_template_size(), toff, tend, caniv[i].template_size())) continue;
			const string& nqstr = op.nqaln();
			const string& ntstr = op.ntaln();
			++num_added;
			used_ids.insert(qid);
			add_one_align(nqstr, ntstr, toff, cov_cnts);
		}
    }
	
	vector<char> cns_seq(op.target_sequence());
    bool has_been_cns = 0;
	int i = 0;
	while (i < template_size) {
		while (i < template_size && cov_cnts[i].cov() < min_cov) ++i;
		int j = i + 1;
		while (j < template_size && cov_cnts[j].cov() >= min_cov) ++j;
		if (j >= template_size) break;
		process_one_frag(cov_cnts, template_id, template_size, i, j, cns_seq);
        has_been_cns = 1;
		i = j;
	}
    if (!has_been_cns) return;
	
	CnsSeq cns_result;
	cns_result.id = extract_read_id_from_ontcns_hdr(op.seq_header(template_id));
	cns_result.left = 0;
	cns_result.right = template_size;
	cns_result.org_seq_size = template_size;
	
	for (size_t k = 0; k != cns_seq.size(); ++k) {
		char c = cns_seq[k];
		r_assert(c >= 0 && c < 5);
		if (c == 4) continue;
		cns_result.seq += "ACGT"[c];
	}
	
	out << cns_result;
	corrected_read_ids.insert(template_id);
}
