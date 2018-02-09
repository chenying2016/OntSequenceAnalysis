#include "error_estimate.h"

#include "../common/oc_assert.h"
#include "../klib/kstring.h"

static BOOL
is_good_overlap(int qoff,
				int qend,
				int qsize,
				int soff,
				int send,
				int ssize)
{
	int qlh = qoff;
	int qrh = qsize - qend;
	int slh = soff;
	int srh = ssize - send;
	const int M = 200;
	
	BOOL r =  (qlh <= M && qrh <= M)
			  ||
			  (slh <= M && srh <= M)
			  ||
			  (qrh <= M && slh <= M)
			  ||
			  (srh <= M && qlh <= M);
	
	return r;
}

static BOOL
estimate_ident_lower_bound(double ident[],
						   const int n,
						   double* avg_error,
						   double* error_cutoff)
{
	if (n < 5) return FALSE;
	
	double sum = 0.0;
	for (int i = 0; i < n; ++i) sum += ident[i];
	double avg = sum / n;
	double se = 0.0;
	for (int i = 0; i < n; ++i) se += (avg - ident[i]) * (avg - ident[i]);
	se /= n;
	se = sqrt(se);
	double s = avg - se * 5;
	
	if (avg_error) *avg_error = avg;
	if (error_cutoff) *error_cutoff = s;
	return TRUE;
}

static void
get_idents(OverlapsPool* op,
		   double* ident,
		   const int NIdent,
		   int* n_ident,
		   const int tsize)
{
	int n = 0;
	for (size_t i = 0; i != oc_op_num_align(*op); ++i) {
		OverlapIndex* oip = oc_op_oip(*op, i);
		if (is_good_overlap(oip->qoff, oip->qend, oip->qsize, oip->toff, oip->tend, tsize)) {
			ident[n++] = oip->ident_perc;
			if (n == NIdent) break;
		}
	}
	if (n < NIdent) {
		n = 0;
		for (size_t i = 0; i != oc_op_num_align(*op); ++i) {
			OverlapIndex* oip = oc_op_oip(*op, i);
			BOOL r = (oip->qend - oip->qoff >= oip->qsize * 0.6) || (oip->tend - oip->toff >= tsize * 0.6);
			if (r) {
				ident[n++] = oip->ident_perc;
				if (n == NIdent) break;
			}
		}
	}
	*n_ident = n;
}

BOOL
get_good_overlaps(GappedCandidate* candidates,
				  size_t can_sid,
				  size_t can_eid,
				  size_t* last_extended_can_id,
				  OcAlignData* align_data,
				  OverlapsPool* op,
				  kstring_t* query,
				  kstring_t* target,
				  PackedDB* reads,
				  ReadIdPool* extended_reads,
				  int min_align_size,
				  double* ident_cutoff)
{
	oc_op_clear(*op);
	clear_ReadIdPool(extended_reads);
	const int NIdent = 10;
	int n_ident = 0;
	double ident[NIdent];
	size_t i;
	for (i = can_sid; i < can_eid && i < can_sid + 50; ++i) {
		GappedCandidate* can = candidates + i;
		//OC_LOG("extending ");
		//DUMP_GAPPED_CANDIDATE(fprintf, stdout, *can);
		if (read_id_exists(extended_reads, can->qid)) continue;
		oc_assert(can->sdir == FWD);
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
		op_add_align(align_data, can->qid, can->qdir, can->qsize, i, op);
		add_read_id(extended_reads, can->qid);
		r = is_good_overlap(oca_query_start(*align_data), 
							oca_query_end(*align_data), 
							can->qsize, 
							oca_target_start(*align_data), 
							oca_target_end(*align_data),
							can->ssize);
		if (r) {
			ident[n_ident++] = oca_ident_perc(*align_data);
			if (n_ident == NIdent) break;
		}
	}
	*last_extended_can_id = i - 1;
	
	if (n_ident < NIdent) get_idents(op, ident, NIdent, &n_ident, candidates[can_sid].ssize);
	return estimate_ident_lower_bound(ident, n_ident, NULL, ident_cutoff);
}
