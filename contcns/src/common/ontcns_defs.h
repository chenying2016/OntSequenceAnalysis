#ifndef ONTCNS_DEFS_H
#define ONTCNS_DEFS_H

/// data type

#include <stdint.h>

typedef int8_t      i8;
typedef uint8_t     u8;
typedef int16_t     i16;
typedef uint16_t    u16;
typedef int32_t     i32;
typedef uint32_t    u32;
typedef int64_t     i64;
typedef uint64_t    u64;

typedef u64         idx;
typedef int8_t      BOOL;

#define IDX_MAX     ((idx)-1)
#define TRUE        1
#define FALSE       0

#define U64_ONE		((u64)1)

inline char
DecodeDNA(const int i)
{
    return "ACGT-"[i];
}

#define FWD 0
#define REV 1
inline int
ReverseStrand(const int s)
{
    return 1 - s;
}

#define GAP_CODE 4
#define GAP_CHAR ('-')

#define TECH_PACBIO 0
#define TECH_NANOPORE 1

#include "kvec.h"

typedef kvec_t(u64) vec_u64;
typedef kvec_t(int) vec_int;
typedef kvec_t(idx) vec_idx;
typedef kvec_t(size_t) vec_size_type;

typedef struct
{
	int first;
	int second;
}IntPair;

typedef kvec_t(IntPair) vec_intpair;
#define IntPair_LT(a, b) (((a).first < (b).first) || ((a).first == (b).first && (a).second < (b).second))
void ks_introsort_IntPair(size_t n, IntPair* a);

#define ARG_PARSE_SUCCESS TRUE
#define ARG_PARSE_FAIL FALSE

#endif // ONTCNS_DEFS_H
