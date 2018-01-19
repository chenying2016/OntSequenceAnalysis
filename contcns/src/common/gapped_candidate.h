#ifndef GAPPED_CANDIDATE_H
#define GAPPED_CANDIDATE_H

#include "ontcns_aux.h"
#include "ontcns_defs.h"
#include "kstring.h"
#include "kvec.h"

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
void ks_introsort_GappedCandidate(size_t n, GappedCandidate* a);

typedef kvec_t(GappedCandidate) vec_can;

void
kstring_output_GappedCandidate(KstringOutputFunc output_func, kstring_t* stream, GappedCandidate* can);

void
buffer_output_GappedCandidate(BufferOutputFunc output_func, char* stream, GappedCandidate* can);

void
stream_output_GappedCandidate(StreamOutputFunc output_func, FILE* stream, GappedCandidate* can);

#endif // GAPPED_CANDIDATE_H
