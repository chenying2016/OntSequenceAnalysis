#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <zlib.h>

#include "kseq.h"
#include "kstring.h"

KSEQ_DECLARE(gzFile)

static inline kstring_t*
kseq_sequence(kseq_t* seq) 
{
	return &seq->seq;
}

static inline kstring_t*
kseq_header(kseq_t* seq)
{
	return &seq->name;
}

static inline size_t
kseq_size(kseq_t* seq)
{
	return KSTR_SIZE(seq->seq);
}

static inline char*
kseq_raw_sequence(kseq_t* seq)
{
	return KSTR_DATA(seq->seq);
}

#endif // SEQUENCE_H
