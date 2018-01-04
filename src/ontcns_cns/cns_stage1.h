#ifndef CNS_STAGE1_H
#define CNS_STAGE1_H

#include "consensus_aux.h"

inline bool
check_overlap_relations(int qoff, int qend, int qsize,
						int toff, int tend, int tsize)
{
    /*
	const int es = 200;
	if (qoff <= es && qsize - qend <= es) return 1;
	if (toff <= es && tsize - tend <= es) return 1;
	if (qoff <= es && tsize - tend <= es) return 1;
	if (toff <= es && qsize - qend <= es) return 1;
	return 0;
    */

    const double R = 0.50;
    bool r = (qend - qoff >= qsize * R) || (tend - toff >= tsize * R);
    if (r) return 1;
    const int A = 2000;
    r = (qend - qoff >= A) || (tend - toff >= A);
    return r;
}

void
cns_stage1(ConsensusData* cns_data,
			const int template_id,
			const int cani_sid,
			const int cani_eid,
			const double min_mapping_ratio);

#endif // CNS_STAGE1_H
