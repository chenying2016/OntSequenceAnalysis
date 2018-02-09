#include "consensus_aux.h"

CnsData*
new_CnsData(PackedDB* reads,
			vec_can* candidates,
			size_t* next_can_id,
			pthread_mutex_t* can_id_lock,
			CnsOptions* options,
			FILE* cns_out,
			pthread_mutex_t* cns_out_lock,
			FILE* raw_out,
			pthread_mutex_t* raw_out_lock,
		    ReadIdPool* corrected_read_ids)
{
	CnsData* data = (CnsData*)malloc(sizeof(CnsData));
	kstr_init(data->query);
	kstr_init(data->target);
	kstr_init(data->qabuf);
	kstr_init(data->tabuf);
	init_CnsSeq(data->cns_seq);
	data->reads = reads;
	data->candidates = candidates;
	data->next_can_id = next_can_id;
	data->can_id_lock = can_id_lock;
	data->options = options;
	data->cns_out = new_ResultsWriter(cns_out, cns_out_lock, FALSE, 0, 64 * (1<<20));
	data->raw_out = new_ResultsWriter(raw_out, raw_out_lock, FALSE, 0, 64 * (1<<20));
	data->cns_data = new_CbCnsData();
	data->align_data = new_OcAlignData(options->error);
	kv_init(data->cov_stats);
	data->op = new_OverlapsPool();
	data->extended_read_ids = new_ReadIdPool(NULL);
	data->corrected_read_ids = corrected_read_ids;
	return data; 
}

CnsData*
free_CnsData(CnsData* data)
{
	free_kstring(data->query);
	free_kstring(data->target);
	free_kstring(data->qabuf);
	free_kstring(data->tabuf);
	free_CnsSeq(data->cns_seq);
	free_ResultsWriter(data->cns_out);
	free_ResultsWriter(data->raw_out);
	free_CbCnsData(data->cns_data);
	free_OcAlignData(data->align_data);
	kv_destroy(data->cov_stats);
	free_OverlapsPool(data->op);
	free_ReadIdPool(data->extended_read_ids);
	free(data);
	return 0;
}

BOOL
extract_candidate_range(CnsData* cns_data, size_t* can_sid, size_t* can_eid)
{
	size_t sid, eid;
	int template_id;
	pthread_mutex_lock(cns_data->can_id_lock);
	sid = *cns_data->next_can_id;
	template_id = kv_A(*cns_data->candidates, sid).sid;
	for (eid = sid + 1; eid < kv_size(*cns_data->candidates); ++eid) {
		if (kv_A(*cns_data->candidates, eid).sid != template_id) break;
	}
	*cns_data->next_can_id = eid;
	pthread_mutex_unlock(cns_data->can_id_lock);
	
	*can_sid = sid;
	*can_eid = eid;
	return sid < kv_size(*cns_data->candidates);
}
