#include "overlaps_pool.h"

#include "../common/ontcns_aux.h"

void
op_add_align(OcAlignData* align_data,
			 const int qid,
			 const int qdir,
			 const int qsize,
			 size_t can_id,
			 OverlapsPool* op)
{
	OverlapIndex oi;
	oi.qid = qid;
	oi.qdir = qdir;
	oi.qoff = oca_query_start(*align_data);
	oi.qend = oca_query_end(*align_data);
	oi.qsize = qsize;
	oi.toff = oca_target_start(*align_data);
	oi.tend = oca_target_end(*align_data);
	oi.ident_perc = oca_ident_perc(*align_data);
	oi.can_id = can_id;
	oi.align_size = kstr_size(*oca_query_mapped_string(*align_data));
	oi.qalign_start = kstr_size(op->align_string);
	oi.talign_start = oi.qalign_start + oi.align_size;
	kv_push(OverlapIndex, op->oiv, oi);
	
	kstring_t* query_align = oca_query_mapped_string(*align_data);
	kputsn(kstr_data(*query_align), kstr_size(*query_align), &op->align_string);
	kstring_t* target_align = oca_target_mapped_string(*align_data);
	kputsn(kstr_data(*target_align), kstr_size(*target_align), &op->align_string);
}

OverlapsPool*
new_OverlapsPool()
{
	OverlapsPool* op = (OverlapsPool*)malloc(sizeof(OverlapsPool));
	kv_init(op->oiv);
	kstr_init(op->align_string);
	return op;
}

OverlapsPool*
free_OverlapsPool(OverlapsPool* op)
{
	kv_destroy(op->oiv);
	free_kstring(op->align_string);
	free(op);
	return 0;
}
