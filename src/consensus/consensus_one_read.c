#include "consensus_one_read.h"

#include "../common/oc_assert.h"
#include "../common/record_writer.h"
#include "../klib/kvec.h"
#include "error_estimate.h"

RESULT_WRITER_IMPL(CnsSeq, CNS_SEQ, static)

static void
get_raw_intvs(int read_size, vec_intpair* cns_intvs, vec_intpair* raw_intvs)
{
	BOOL first = TRUE;
	BOOL last_interval_open = TRUE;
	int start_offset, end_offset, filter_start, filter_end;
	int left = 0, right;
	IntPair intv;
	
	start_offset = 0;
	end_offset = read_size - 1;
	for (size_t i = 0; i != kv_size(*cns_intvs); ++i) {
		filter_start = start_offset + kv_A(*cns_intvs, i).first;
		filter_end = start_offset + kv_A(*cns_intvs, i).second;
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
		if (right - left + 1 >= 1000) kv_push(IntPair, *raw_intvs, intv);
		
		if (filter_end >= end_offset) {
			last_interval_open = 0;
			break;
		} else {
			left = filter_end + 1;
		}
	}
	
	if (last_interval_open) {
		right = end_offset;
		if (right - left + 1 >= 1000) {
			intv.first = left;
			intv.second = right;
			kv_push(IntPair, *raw_intvs, intv);
		}
	}
}

static void
consensus_one_read_broken(CbCnsData* cbcns_data,
						 kstring_t* target,
						 const int min_cov,
						 const int min_size,
						 const int template_id,
						 const int template_size,
						 CnsSeq* cns_seq,
						 ResultsWriter* cns_out,
						 ResultsWriter* raw_out)
{
	new_kvec(vec_intpair, cns_intvs);
	new_kvec(vec_intpair, raw_intvs);
	consensus_broken(cbcns_data,
					  min_cov,
					  min_size,
					  template_id,
					  template_size,
					  &cns_intvs,
					  cns_seq,
					  cns_out);
	get_raw_intvs(template_size, &cns_intvs, &raw_intvs);
	for (size_t i = 0; i != kv_size(raw_intvs); ++i) {
		int from = kv_A(raw_intvs, i).first;
		int to = kv_A(raw_intvs, i).second + 1;
		cns_seq->id = template_id;
		cns_seq->left = from;
		cns_seq->right = to;
		cns_seq->org_seq_size = template_size;
		kstr_clear(cns_seq->cns_seq);
		for (int k = from; k < to; ++k) {
			char c = DecodeDNA( kstr_A(*target, k) );
			kputc(c, &cns_seq->cns_seq);
		}
		RW_DUMP_ONE_DATA(CNS_SEQ, raw_out, cns_seq);
	}
	free_kvec(cns_intvs);
	free_kvec(raw_intvs);
}

static void
consensus_one_read_unbroken(CbCnsData* cbcns_data,
							  kstring_t* target,
							  const int min_cov,
							  const int min_size,
							  const int template_id,
							  const int template_size,
							  CnsSeq* cns_seq,
							  ResultsWriter* cns_out,
							  ResultsWriter* raw_out)
{
	consensus_unbroken(cbcns_data, min_cov, min_size, kstr_str(*target), template_size, &cns_seq->cns_seq);
	cns_seq->id = template_id;
	cns_seq->left = 0;
	cns_seq->right = template_size;
	cns_seq->org_seq_size = template_size;
	RW_DUMP_ONE_DATA(CNS_SEQ, cns_out, cns_seq);
}

void
update_cov_stats(int* cov_stats, int toff, int tend)
{
	for (int i = toff; i < tend; ++i) ++cov_stats[i];
}

BOOL
region_coverage_is_full(const int* cov_stats, int toff, int tend, int max_cov)
{
	int n = 0;
	for (int i = toff; i < tend; ++i) if (cov_stats[i] >= max_cov) ++n;
	return n + 200 >= tend - toff;
}

void
add_extended_overlaps(OverlapsPool* op,
					  vec_int* cov_stats,
					  ReadIdPool* extended_read_ids,
					  const double weight,
					  const double ident_cutoff,
					  kstring_t* qabuf,
					  kstring_t* tabuf,
					  CbCnsData* cbcns_data,
					  kstring_t* target)
{
	for (size_t i = 0; i != oc_op_num_align(*op); ++i) {
		OverlapIndex* oip = oc_op_oip(*op, i);
		if (oip->ident_perc < ident_cutoff) continue;
		add_one_align(cbcns_data, oc_op_qstr(*op, i), oc_op_tstr(*op, i), oip->align_size, target, oip->toff, oip->tend, weight);
		update_cov_stats(kv_data(*cov_stats), oip->toff, oip->tend);
		add_read_id(extended_read_ids, oip->qid);
	}
}

void
consensus_one_read(CnsData* cns_data, size_t can_sid, size_t can_eid)
{
	if (cns_data->options->high_accuracy_coverage_cutoff > can_eid - can_sid) return;
	
	PackedDB* reads = cns_data->reads;
	GappedCandidate* candidates = kv_data(*cns_data->candidates);
	CnsOptions* options = cns_data->options;
	ResultsWriter* cns_out = cns_data->cns_out;
	ResultsWriter* raw_out = cns_data->raw_out;
	ReadIdPool* corrected_read_ids = cns_data->corrected_read_ids;
	ReadIdPool* extended_read_ids = cns_data->extended_read_ids;
	CbCnsData* cbcns_data = cns_data->cns_data;
	OcAlignData* align_data = cns_data->align_data;
	vec_int* cov_stats = &cns_data->cov_stats;
	OverlapsPool* op = cns_data->op;
	const int min_align_size = options->align_size_cutoff;
	const int min_size = options->min_size;
	const int template_id = candidates[can_sid].sid;
	const int template_size = candidates[can_sid].ssize;
	kstring_t* query = &cns_data->query;
	kstring_t* target = &cns_data->target;
	kstring_t* qabuf = &cns_data->qabuf;
	kstring_t* tabuf = &cns_data->tabuf;
	CnsSeq* cns_seq = &cns_data->cns_seq;
	
	//if (template_id != 226) return;
	//OC_LOG("consensus can %lu --- %lu, read_id = %d, read_size = %d", can_sid, can_eid, template_id, template_size);
	
	pdb_extract_sequence(reads, template_id, FWD, target);
	ks_introsort_GappedCandidate_ScoreGT(can_eid - can_sid, candidates + can_sid);
	clear_CbCnsData(cbcns_data);
	cbcns_data->template_size = template_size;
	kv_resize(int, *cov_stats, template_size);
	kv_zero(int, *cov_stats);
	
	size_t last_extended_can_id;
	double ident_cutoff;
	const double cns_weight = 1.0;
	if (!get_good_overlaps(candidates,
						   can_sid,
						   can_eid,
						   &last_extended_can_id,
						   align_data,
						   op,
						   query,
						   target,
						   reads,
						   extended_read_ids,
						   min_align_size,
						   &ident_cutoff)) {
		return;
	}
	//OC_LOG("ident cutoff: %f", ident_cutoff);
	
	const int max_cov = (ident_cutoff >= 0.70) ? options->high_accuracy_coverage_cutoff : options->low_accuracy_coverage_cutoff;
	add_extended_overlaps(op, cov_stats, extended_read_ids, cns_weight, ident_cutoff, qabuf, tabuf, cbcns_data, target);
	
	size_t next_can_id = last_extended_can_id + 1;
	while (1) {
		if (next_can_id >= can_eid) break;
		if (region_coverage_is_full(kv_data(*cov_stats), 0, template_size, max_cov)) break;
		size_t from = next_can_id;
		size_t to = from + 50; to = OC_MIN(to, can_eid);
		next_can_id = to;
		for (size_t i = from; i < to; ++i) {
			GappedCandidate* can = candidates + i;
			oc_assert(can->sdir == FWD);
			oc_assert(can->sid == template_id);
			if (read_id_exists(extended_read_ids, can->qid)) continue;
			if (region_coverage_is_full(kv_data(*cov_stats), can->sbeg, can->send, max_cov)) continue;
			pdb_extract_sequence(reads, can->qid, can->qdir, query);
			align_data->qid = can->qid;
			align_data->tid = can->sid;
			BOOL r = onc_align(kstr_str(*query),
							   can->qoff,
							   can->qsize,
							   kstr_str(*target),
							   can->soff,
							   can->ssize,
							   align_data,
							   min_align_size);
			if (!r) continue;
			if (oca_ident_perc(*align_data) < ident_cutoff) continue;
			kstring_t* qalign = oca_query_mapped_string(*align_data);
			kstring_t* talign = oca_target_mapped_string(*align_data);
			add_one_align(cbcns_data, kstr_data(*qalign), kstr_data(*talign), kstr_size(*qalign), target, oca_target_start(*align_data), oca_target_end(*align_data), cns_weight);
			update_cov_stats(kv_data(*cov_stats), oca_target_start(*align_data), oca_target_end(*align_data));
			add_read_id(extended_read_ids, can->qid);
		}
	}
	
	if (options->full_consensus) {
		consensus_one_read_unbroken(cbcns_data, target, max_cov, min_size, template_id, template_size, cns_seq, cns_out, raw_out);
	} else {
		consensus_one_read_broken(cbcns_data, target, max_cov, min_size, template_id, template_size, cns_seq, cns_out, raw_out);
	}
	
	add_read_id(corrected_read_ids, template_id);
}
