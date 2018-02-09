#ifndef CONSENSUS_AUX_H
#define CONSENSUS_AUX_H

#include "../common/packed_db.h"
#include "../common/gapped_candidate.h"
#include "../common/record_writer.h"
#include "../gapped_align/oc_aligner.h"
#include "../tasc/cbcns.h"
#include "cns_options.h"
#include "overlaps_pool.h"
#include "read_id_pool.h"

typedef struct {
	kstring_t			query;
	kstring_t			target;
	kstring_t			qabuf;
	kstring_t			tabuf;
	CnsSeq				cns_seq;
	PackedDB* 			reads;
	vec_can*			candidates;
	size_t*				next_can_id;
	pthread_mutex_t*	can_id_lock;
	CnsOptions*			options;
	ResultsWriter*		cns_out;
	ResultsWriter*		raw_out;
	CbCnsData* 			cns_data;
	OcAlignData* 		align_data;
	vec_int				cov_stats;
	OverlapsPool* 		op;
	ReadIdPool*			extended_read_ids;
	ReadIdPool*			corrected_read_ids;
} CnsData;

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
		    ReadIdPool* corrected_read_ids);

CnsData*
free_CnsData(CnsData* data);

BOOL
extract_candidate_range(CnsData* cns_data, size_t* can_sid, size_t* can_eid);

#endif // CONSENSUS_AUX_H
