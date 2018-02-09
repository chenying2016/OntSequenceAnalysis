#ifndef GAPPED_CANDIDATE_H
#define GAPPED_CANDIDATE_H

#include "ontcns_aux.h"
#include "ontcns_defs.h"
#include "../klib/kstring.h"
#include "../klib/kvec.h"

typedef struct 
{
	int qid;
	int sid;
	int qdir;
	int sdir;
	int score;
	idx qbeg, qend, qsize;
	idx sbeg, send, ssize;
	idx qoff, soff;
} GappedCandidate;

#define GappedCandidate_ScoreGT(a, b) ((a).score > (b).score)
void ks_introsort_GappedCandidate_ScoreGT(size_t n, GappedCandidate* a);

#define GappedCandidate_SidLT(a, b) ((a).sid < (b).sid)
void ks_introsort_GappedCandidate_SidLT(size_t n, GappedCandidate* a);

typedef kvec_t(GappedCandidate) vec_can;

#define DUMP_GAPPED_CANDIDATE(output_func, out, can) \
	output_func(out, \
			"%d\t%d\t%d\t%d\t%lu\t%lu\t%lu\t%lu\t%d\t%lu\t%lu\t%lu\t%lu\n", \
			(can).qid, \
			(can).sid, \
			(can).score, \
			(can).qdir, \
			(can).qbeg, \
			(can).qend, \
			(can).qoff, \
			(can).qsize, \
			(can).sdir, \
			(can).sbeg, \
			(can).send, \
			(can).soff, \
			(can).ssize) \

#define LOAD_GAPPED_CANDIDATE(input_func, in, can) \
	input_func(in, "%d\t%d\t%d\t%d\t%lu\t%lu\t%lu\t%lu\t%d\t%lu\t%lu\t%lu\t%lu\n", \
					&(can).qid, \
					&(can).sid, \
					&(can).score, \
					&(can).qdir, \
					&(can).qbeg, \
					&(can).qend, \
					&(can).qoff, \
					&(can).qsize, \
					&(can).sdir, \
					&(can).sbeg, \
					&(can).send, \
					&(can).soff, \
					&(can).ssize)

#endif // GAPPED_CANDIDATE_H
