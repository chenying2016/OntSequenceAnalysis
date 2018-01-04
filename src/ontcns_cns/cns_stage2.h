#ifndef CNS_STAGE2_H
#define CNS_STAGE2_H

#include "consensus_aux.h"

void
cns_stage2(ConsensusData* cns_data,
			const int template_id,
			const int cani_sid,
			const int cani_eid,
			const double min_mapping_ratio);

void
cns_stage2(ConsensusData* cns_data,
		   size_t can_sid,
		   size_t can_eid);

#endif // CNS_STAGE2_H
