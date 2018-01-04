#include "cns_stage2.h"

#include "../common/cns_seq.h"
#include "../common/smart_assert.h"
#include "consensus_aux.h"
#include "cns_stage1.h"

#include <algorithm>
#include <set>
#include <string>

using namespace std;

static const int MaxCov = 40;
static const int MaxUnCnsSize = 200;

double 
calc_cns_weight(const double error)
{
	double w = (1.0 - error) * (1.0 - error) + error * error / 3.0;
	return w;
}

bool
cns_full(CntBaseConsensus& cbcns,
		 const int min_cov,
		 const int min_size,
		 const int template_id,
		 const vector<char>& raw_read,
		 OutputStream& out)
{
	CnsSeq cns_seq;
	cns_seq.id = template_id;
	cns_seq.left = 0;
	cns_seq.right = raw_read.size();
	cns_seq.org_seq_size = raw_read.size();
	if (!cbcns.consensus(min_cov, min_size, raw_read, cns_seq.seq)) return 0;
	r_assert(cns_seq.seq.size());
	out << cns_seq;
	return 1;
}

void
get_raw_intvs(int read_size,
			  vector<pair<int, int> >& cns_intvs,
			  vector<pair<int, int> >& raw_intvs)
{
	bool first = 1;
	bool last_interval_open = 1;
	int start_offset, end_offset, filter_start, filter_end;
	int left = 0, right;
	pair<int, int> intv;
	
	start_offset = 0;
	end_offset = read_size - 1;
	for (size_t i = 0; i != cns_intvs.size(); ++i) {
		filter_start = start_offset + cns_intvs[i].first;
		filter_end = start_offset + cns_intvs[i].second;
		if (first) {
			last_interval_open = 0;
			first = 0;
			if (filter_start > start_offset) {
				left = start_offset;
			} else {
				left = filter_end + 1;
				continue;
			}
		}
		
		right = filter_start - 1;
		intv.first = left;
		intv.second = right;
		if (right - left + 1 >= 1000) raw_intvs.push_back(make_pair(left, right));
		
		if (filter_end >= end_offset) {
			last_interval_open = 0;
			break;
		} else {
			left = filter_end + 1;
		}
	}
	
	if (last_interval_open) {
		right = end_offset;
		if (right - left + 1 >= 1000) raw_intvs.push_back(make_pair(left, right));
	}
}

bool
cns_nonfull(CntBaseConsensus& cbcns,
			const vector<char>& target,
			const int min_cov,
			const int min_size,
			const int template_id,
			const int template_size,
			OutputStream& out,
		    OutputStream& raw_out)
{
	vector<pair<int, int> > cns_intvs, raw_intvs;
	if (!cbcns.consensus(min_cov, min_size, template_id, template_size, cns_intvs, out)) return 0;
	get_raw_intvs(template_size, cns_intvs, raw_intvs);
	CnsSeq cns_seq;
	for (size_t i = 0; i != raw_intvs.size(); ++i) {
		int from = raw_intvs[i].first;
		int to = raw_intvs[i].second + 1;
		
		cns_seq.id = template_id;
		cns_seq.left = from;
		cns_seq.right = to;
		cns_seq.org_seq_size = template_size;
		cns_seq.seq.clear();
		for (int k = from; k < to; ++k) {
			cns_seq.seq += DNADecoder::decode( target[k] ); 
		}
		raw_out << cns_seq;
	}
	return 1;
}

bool
qualify_can_cov(vector<int>& cov_stats, GappedCandidate& c)
{
	int n = 0;
	for (int i = c.sbeg; i < c.send; ++i) {
		if (cov_stats[i] < MaxCov) ++n;
	}
	return n >= MaxUnCnsSize;
}

bool
check_cov_stats(vector<int>& cov_stats) 
{
	int n = 0;
	for (size_t i = 0; i != cov_stats.size(); ++i) {
		if (cov_stats[i] < MaxCov) ++n;
	}
	return n >= MaxUnCnsSize;
}

bool
find_paired_alignment(multimap<int, size_t>& qids,
					  const int qid,
					  const int qdir,
					  const int qoff,
					  const int qend,
					  const int qsize,
					  const int toff,
					  const int tend,
					  const int tsize,
					  const size_t can_id,
					  const int min_align_size,
					  vector<char>& query,
					  vector<char>& target,
					  vector<GappedCandidate>& candidates,
					  OntCnsAligner& aligner)
{
	const int MinHang = 500;
	bool lhang = (qoff >= MinHang) && (toff >= MinHang);
	bool rhang = (qsize - qend >= MinHang) && (tsize - tend >= MinHang);
	if ((!lhang) && (!rhang)) return 0;
	const size_t kInvalidSize = 100000000;
	size_t lcid = kInvalidSize;
	size_t rcid = kInvalidSize;
	
	pair< multimap<int, size_t>::iterator, multimap<int, size_t>::iterator > range = qids.equal_range(qid);
	for (multimap<int, size_t>::iterator iter = range.first; iter != range.second; ++iter) {
		if (iter->second == can_id) continue;
		r_assert(iter->first == qid);
		GappedCandidate& can = candidates[iter->second];
		if (can.qdir != qdir) continue;
		if (lhang) {
			if (can.qoff < qoff && can.soff < toff) {
				lcid = iter->second;
				break;
			}
		} 
		if (rhang) {
			if (can.qoff > qend && can.soff > tend) {
				rcid = iter->second;
				break;
			}
		}
	}
	
	size_t cid = (lcid == kInvalidSize) ? rcid : lcid;
	if (cid == kInvalidSize) return 0;
	GappedCandidate& can = candidates[cid];
	aligner.set_ids(can.qid, can.sid);
	bool r = aligner.go(query.data(), can.qoff, can.qsize, target.data(), can.soff, can.ssize, min_align_size);
	return r;
}

void
update_cov_stats(vector<int>& cov_stats, int soff, int send)
{
	for (int i = soff; i < send; ++i) ++cov_stats[i];
}

bool
is_good_overlap(int qoff,
				int qend,
				int qsize,
				int soff,
				int send,
				int ssize)
{
	int qlh = qoff;
	int qrh = qsize - qend;
	int tlh = soff;
	int trh = ssize - send;
	
	const int M = 200;
	if (qlh <= M && qrh <= M) return 1;
	if (tlh <= M && trh <= M) return 1;
	
	if (qrh <= M && tlh <= M) return 1;
	if (trh <= M && qlh <= M) return 1;
	
	return 0;
}

bool
calc_ident_lower_bound(double ident[], 
					   const int n, 
					   double& avg_error,
					   double& error_cutoff)
{
	if (n < 5) return 0;
	
	double sum = 0.0;
	for (int i = 0; i < n; ++i) sum += ident[i];
	double avg = sum / n;
	double se = 0;
	for (int i = 0; i < n; ++i) se += (avg - ident[i]) * (avg - ident[i]);
	se /= (n);
	se = sqrt(se);
	double s =  avg - se * 5;
	avg_error = avg;
	error_cutoff = s;
	return 1;
}

void
get_idents(OverlapInfoV& oiv,
		   double* ident,
		   const int NIdent,
		   int& n_ident)
{
	n_ident = 0;
	for (size_t i = 0; i != oiv.oiv.size(); ++i) {
		OverlapInfo& oi = oiv.oiv[i];
		if (is_good_overlap(oi.qoff, oi.qend, oi.qsize, oi.toff, oi.tend, oi.tsize)) {
			ident[n_ident++] = oi.ident_perc;
			if (n_ident == NIdent) break;
		}
	}
	if (n_ident < NIdent) {
		n_ident = 0;
		for (size_t i = 0; i != oiv.oiv.size(); ++i) {
			OverlapInfo& oi = oiv.oiv[i];
			bool r = (oi.qend - oi.qoff >= oi.qsize * 0.6) || (oi.tend - oi.toff >= oi.tsize * 0.6);
			if (r) {
				ident[n_ident++] = oi.ident_perc;
				++n_ident;
			}
			if (n_ident == NIdent) break;
		}
	}
}

void
get_good_overlaps(vector<GappedCandidate>& candidates,
				  size_t can_sid,
				  size_t can_eid,
				  size_t& last_extended_can_id,
				  OntCnsAligner& aligner,
				  OverlapInfoV& oiv_pb,
				  OverlapInfoV& oiv_ont,
				  vector<char>& query,
				  vector<char>& target,
				  PackedDB& reads,
				  int min_align_size,
				  double& pb_ident_cutoff,
				  double& pb_weight,
				  bool& pb_data_effective,
				  double& ont_ident_cutoff,
				  double& ont_weight,
				  bool& ont_data_effective)
{
	oiv_pb.clear();
	oiv_ont.clear();
	size_t i;
	int num_pb_can = 0;
	int num_ont_can = 0;
	set<int> used_ids;
	for (i = can_sid; i < can_eid && i < can_sid + 100; ++i) {
		GappedCandidate& can = candidates[i];
		if (used_ids.find(can.qid) != used_ids.end()) continue;
		query.resize(can.qsize);
		reads.get_sequence(can.qid, can.qdir == FWD, query.data());
		int tech = reads.seq_tech(can.qid);
		aligner.set_ids(can.qid, can.sid);
		if (!aligner.go(query.data(), can.qoff, can.qsize, target.data(), can.soff, can.ssize, min_align_size)) continue;
		if (tech == TECH_PACBIO) {
			++num_pb_can;
			oiv_pb.add(aligner, can.qid, can.qdir, can.qsize, can.ssize, i);
		} else if (tech == TECH_NANOPORE) {
			++num_ont_can;
			oiv_ont.add(aligner, can.qid, can.qdir, can.qsize, can.ssize, i);
		} else {
			mc_error << "invalid Sequencing Platform: " << tech << eolog;
		}
		used_ids.insert(can.qid);
	}
	last_extended_can_id = i - 1;
	
	const int NIdent = 10;
	int n_ident = 0;
	double ident[NIdent];
	
	pb_data_effective = false;
	if (num_pb_can) {
		get_idents(oiv_pb, ident, NIdent, n_ident);
		if (calc_ident_lower_bound(ident, n_ident, pb_weight, pb_ident_cutoff)) {
			pb_weight = calc_cns_weight(pb_weight);
			pb_data_effective = true;
		}
	}
	
	ont_data_effective = false;
	if (num_ont_can) {
		get_idents(oiv_ont, ident, NIdent, n_ident);
		if (calc_ident_lower_bound(ident, n_ident, ont_weight, ont_ident_cutoff)) {
			ont_weight = calc_cns_weight(ont_weight);
			ont_data_effective = true;
		}
	}
}

void
add_preextended_ovlps(vector<GappedCandidate>& candidates,
					  OverlapInfoV& oiv,
					  vector<int>& cov_stats,
					  CntBaseConsensus& cbcns,
					  vector<char>& query,
					  vector<char>& target,
					  PackedDB& reads,
					  OntCnsAligner& aligner,
					  set<int>& used_ids,
					  multimap<int, size_t>& qids,
					  const double weight,
					  const double ident_cutoff,
					  const int min_align_size,
					  const bool need_normalise_gaps)
{
	string qalign, talign;
	for (size_t i = 0; i != oiv.oiv.size(); ++i) {
		OverlapInfo& oi = oiv.oiv[i];
		if (oi.ident_perc < ident_cutoff) continue;
		qalign.assign(oiv.align_string.data() + oi.qalign_start, oi.align_size);
		talign.assign(oiv.align_string.data() + oi.talign_start, oi.align_size);
		cbcns.add_one_align(qalign, talign, oi.toff, weight, need_normalise_gaps);
		update_cov_stats(cov_stats, oi.toff, oi.tend);
		query.resize(oi.qsize);
		reads.get_sequence(oi.qid, oi.qdir == FWD, query.data());
		bool r = find_paired_alignment(qids,
									   oi.qid,
									   oi.qdir,
									   oi.qoff,
									   oi.qend,
									   oi.qsize,
									   oi.toff,
									   oi.tend,
									   oi.tsize,
									   oi.can_id,
									   min_align_size,
									   query,
									   target,
									   candidates,
									   aligner);
		used_ids.insert(oi.qid);
		if (r) {
			cbcns.add_one_align(aligner.query_mapped_string(), 
								aligner.target_mapped_string(),
								aligner.target_start(),
								weight,
							    need_normalise_gaps);
			update_cov_stats(cov_stats, aligner.target_start(), aligner.target_end());
		}
	}
}

void
cns_stage2(ConsensusData* cns_data,
		   size_t can_sid,
		   size_t can_eid)
{
	if (can_eid - can_sid < cns_data->options->min_cov) return;
	
	PackedDB& reads = *cns_data->reads;
	vector<GappedCandidate>& candidates = *cns_data->candidates;
	CnsOptions& options = *cns_data->options;
	OutputStream& out = cns_data->out;
	OutputStream& raw_out = cns_data->raw_out;
	set<int>& corrected_read_ids = cns_data->corrected_read_ids;
	CntBaseConsensus& cbcns = cns_data->cbcns;
	OntCnsAligner& aligner = cns_data->aligner;
	OntCnsAligner& aligner2 = cns_data->aligner2;
	vector<int>& cov_stats = cns_data->cov_stats;
	OverlapInfoV& oiv_pb = cns_data->oiv_pb;
	OverlapInfoV& oiv_ont = cns_data->oiv_ont;
	const int min_align_size = options.min_align_size;
	const int min_cov = options.min_cov;
	const int min_size = options.min_size;
	const int template_id = candidates[can_sid].sid;
	const int template_size = candidates[can_sid].ssize;
	vector<char> target;
	vector<char> query;
	multimap<int, size_t> qids;
	set<int> used_ids;
	bool need_normalise_gaps = options.normalise_gaps;
	
	target.resize(template_size);
	reads.get_sequence(template_id, 1, target.data());
	sort(candidates.begin() + can_sid, candidates.begin() + can_eid, GappedCandidateScore_GT());
	for (size_t i = can_sid; i < can_eid; ++i) {
		r_assert(candidates[i].sid == template_id);
		qids.insert( make_pair(candidates[i].qid, i) );
	}
	cbcns.clear();
	cov_stats.resize(template_size);
	fill(cov_stats.begin(), cov_stats.end(), 0);
	
	size_t last_extended_can_id;
	double pb_ident_cutoff, ont_ident_cutoff;
	double pb_weight, ont_weight;
	bool pb_data_effective, ont_data_effective;
	get_good_overlaps(candidates, 
					  can_sid, 
					  can_eid, 
					  last_extended_can_id, 
					  aligner, 
					  oiv_pb, 
					  oiv_ont, 
					  query, 
					  target, 
					  reads, 
					  min_align_size,
					  pb_ident_cutoff,
					  pb_weight,
					  pb_data_effective,
					  ont_ident_cutoff,
					  ont_weight,
					  ont_data_effective);
	
	if ((!pb_data_effective) && (!ont_data_effective)) return;
	if (!pb_data_effective) {
		pb_ident_cutoff = 2.0;
		pb_weight = 0.0;
		ont_weight = 1.0;
	}
	if (!ont_data_effective) {
		ont_ident_cutoff = 2.0;
		ont_weight = 0.0;
		pb_weight = 1.0;
	}
	
	add_preextended_ovlps(candidates,
						  oiv_pb,
						  cov_stats,
						  cbcns,
						  query,
						  target,
						  reads,
						  aligner,
						  used_ids,
						  qids,
						  pb_weight,
						  pb_ident_cutoff,
						  min_align_size,
						  need_normalise_gaps);
	
	add_preextended_ovlps(candidates,
						  oiv_ont,
						  cov_stats,
						  cbcns,
						  query,
						  target,
						  reads,
						  aligner,
						  used_ids,
						  qids,
						  ont_weight,
						  ont_ident_cutoff,
						  min_align_size,
						  need_normalise_gaps);
	
	{
		for (size_t i = last_extended_can_id + 1; i < can_eid && i < can_sid + 500; ++i) {
			GappedCandidate& can = candidates[i];
			int qid = can.qid;
			if (used_ids.find(qid) != used_ids.end()) continue;
			if (!qualify_can_cov(cov_stats, can)) continue;
			query.resize(can.qsize);
			reads.get_sequence(can.qid, can.qdir == FWD, query.data());
			aligner.set_ids(can.qid, can.sid);
			if (!aligner.go(query.data(), can.qoff, can.qsize, target.data(), can.soff, can.ssize, min_align_size)) continue;
			int tech = reads.seq_tech(can.qid);
			double weight;
			double ident_cutoff;
			if (tech == TECH_PACBIO) {
				weight = pb_weight;
				ident_cutoff = pb_ident_cutoff;
			} else {
				weight = ont_weight;
				ident_cutoff = ont_ident_cutoff;
			}
			if (aligner.calc_ident_perc() < ident_cutoff) continue;
			cbcns.add_one_align(aligner.query_mapped_string(),
								aligner.target_mapped_string(),
								aligner.target_start(),
								weight, 
								need_normalise_gaps);
			update_cov_stats(cov_stats, aligner.target_start(), aligner.target_end());
			used_ids.insert(qid);
			bool r = find_paired_alignment(qids,
										can.qid,
										can.qdir,
										aligner.query_start(),
										aligner.query_end(),
										can.qsize,
										aligner.target_start(),
										aligner.target_end(),
										can.ssize,
										i,
										min_align_size,
										query,
										target,
										candidates,
										aligner2);
			if (r) {
				cbcns.add_one_align(aligner2.query_mapped_string(), 
									aligner2.target_mapped_string(),
									aligner2.target_start(),
									weight,
									need_normalise_gaps);
				update_cov_stats(cov_stats, aligner2.target_start(), aligner2.target_end());
			}
		} // for
	} 
	
	bool has_been_cns = 0;
	const int gid = extract_read_id_from_ontcns_hdr(reads.seq_header(template_id));
	if (cns_data->options->full_cns) {
		has_been_cns = cns_full(cbcns,
								min_cov,
								min_size,
								gid,
								target,
								out);
	} else {
		has_been_cns = cns_nonfull(cbcns,
								   target,
								   min_cov,
								   min_size,
								   gid,
								   template_size,
								   out,
								   raw_out);
	}
	if (has_been_cns) corrected_read_ids.insert(template_id);
}
