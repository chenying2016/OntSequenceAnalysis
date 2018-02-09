#ifndef EDLIB_EX_AUX_H
#define EDLIB_EX_AUX_H

#include <stdint.h>

#include "../common/ontcns_defs.h"
#include "../klib/kstring.h"
#include "../klib/kvec.h"

typedef uint64_t        Word;
#define WORD_SIZE       64
#define WORD_ONE        ((Word)1)
#define HIGH_BIT_MASK   (WORD_ONE << (WORD_SIZE - 1))

#define MaxSeqSize      1024
#define MaxNumBlocks    16
#define AlphabetSize    4

#define EDLIB_OP_SUB 	('M')
#define EDLIB_OP_INS 	('I')
#define EDLIB_OP_DEL 	('D')

/**
 * Describes cigar format.
 * @see http://samtools.github.io/hts-specs/SAMv1.pdf
 * @see http://drive5.com/usearch/manual/cigar.html
 */
typedef enum{
	EDLIB_CIGAR_STANDARD,  //!< Match: 'M', Insertion: 'I', Deletion: 'D', Mismatch: 'M'.
	EDLIB_CIGAR_EXTENDED   //!< Match: '=', Insertion: 'I', Deletion: 'D', Mismatch: 'X'.
}EdlibCigarFormat;

typedef enum 
{
    EDLIB_MODE_NW,
    EDLIB_MODE_SHW,
    EDLIB_MODE_HW
}EdlibAlignMode;

typedef struct {
	Word* Ps;
	Word* Ms;
	int* scores;
	int* first_blocks;
	int* last_blocks;
}EdlibAlignMatrix;

EdlibAlignMatrix*
new_EdlibAlignMatrix();

EdlibAlignMatrix*
free_EdlibAlignMatrix(EdlibAlignMatrix* data);

typedef struct {
	Word P;
	Word M;
	int score;
}EdlibBlock;

typedef struct {
	int 		edit_distance;
	vec_int		end_locations;
	vec_int 	start_locations;
	vec_uchar	alignment;
}EdlibAlignResult;

EdlibAlignResult*
new_EdlibAlignResult();

EdlibAlignResult*
free_EdlibAlignResult(EdlibAlignResult* result);

void
clear_EdlibAlignResult(EdlibAlignResult* result);

typedef struct {
	int num;
	int op;
} EdlibGapAlignOp;

typedef kvec_t(EdlibGapAlignOp)	vec_EdlibGapOp;

typedef struct {
    EdlibAlignMatrix*   align_matrix;
    EdlibAlignResult*   result;
    EdlibBlock*         blocks;
    Word*               peq;
    vec_EdlibGapOp      cigar;
} EdlibAlignData;

EdlibAlignData*
new_EdlibAlignData();

EdlibAlignData*
free_EdlibAlignData(EdlibAlignData* data);

#endif // EDLIB_EX_AUX_H
