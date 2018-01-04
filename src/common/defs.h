//
// Created by sysu on 2/26/17.
//

#ifndef UNTITLED_DEFS_H
#define UNTITLED_DEFS_H

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

#define I8_MIN      ((i8)-128)
#define I16_MIN     ((i16)-32768)
#define I32_MIN     ((i32)-2147483648)
#define I64_MIN     ((i64)-9223372036854775808)

#define I8_MAX      ((i8)127)
#define I16_MAX     ((i16)32767)
#define I32_MAX     ((i32)2147483647)
#define I64_MAX     ((i64)9223372036854775807)

#define U8_MAX      ((u8)255)
#define U16_MAX     ((u16)65535)
#define U32_MAX     ((u32)4294967295)
#define U64_MAX     ((u64)18446744073709551615)

typedef i64         idx;
#define IDX_MIN     I64_MIN
#define IDX_MAX     I64_MAX
typedef u64         uidx;
#define UIDX_MIN    ((u64)0)
#define UIDX_MAX    U64_MAX
#define INVALID_IDX IDX_MAX

#define BEGIN_NAMESPACE_SCOPE {
#define END_NAMESPACE_SCOPE }

/// dna sequence

const u8* get_dna_encode_table();
const char* get_dna_decode_table();
const u8* get_dna_complement_table();

struct DNADecoder
{
static char decode(const int i) {
	return "ACGT-"[i];
}
};

#define FWD 0
#define REV 1
#define REVERSE_STRAND(s) (1 - (s))
#define MAX_SEQ_SIZE (500000)

#define SANITY_CHECK 0
#define SENTINEL (-1)
#define GAP_CODE 4
#define GAP_CHAR ('-')

#define TECH_PACBIO 0
#define TECH_NANOPORE 1

#endif //UNTITLED_DEFS_H
