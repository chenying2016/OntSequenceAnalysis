#include "gapped_candidate.h"

#include "ksort.h"

KSORT_INIT(GappedCandidate, GappedCandidate, GappedCandidate_ScoreGT)

#define CAN_OUTPUT_IMPL \
	output_func(stream, \
				"%d\t%d\t%d\t%d\t%lu\t%lu\t%lu\t%lu\t%d\t%lu\t%lu\t%lu\t%lu\n", \
				can->qid, \
				can->sid, \
				can->score, \
				can->qdir, \
				can->qbeg, \
				can->qend, \
				can->qoff, \
				can->qsize, \
				can->sdir, \
				can->sbeg, \
				can->send, \
				can->soff, \
				can->ssize) \

void
kstring_output_GappedCandidate(KstringOutputFunc output_func, kstring_t* stream, GappedCandidate* can)
{
	CAN_OUTPUT_IMPL;
}

void
buffer_output_GappedCandidate(BufferOutputFunc output_func, char* stream, GappedCandidate* can)
{
	CAN_OUTPUT_IMPL;
}

void
stream_output_GappedCandidate(StreamOutputFunc output_func, FILE* stream, GappedCandidate* can)
{
	CAN_OUTPUT_IMPL;
}
