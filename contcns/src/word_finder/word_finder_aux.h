#ifndef WORD_FINDER_AUX_H
#define WORD_FINDER_AUX_H

#include "../common/kvec.h"
#include "../common/ontcns_defs.h"

#define SM 40
#define SI (SM + 1)

typedef struct
{
    idx qoff;
    idx soff;
} ChainSeed;

typedef kvec_t(ChainSeed) vec_chain_seed;
#define ChainSeedLT(a, b) ((a).soff < (b).soff || ((a).soff == (b).soff && (a).qoff < (b).qoff))

typedef struct
{
	int qno;
	int boff;
} BlockSeed;

typedef struct 
{
	BlockSeed bseeds[SM];
	short score;
	short score2;
} BlockSeeds;

typedef struct
{
	idx start;
	idx cnt;
	idx block_id;
} BlockSeedsIndex;

typedef kvec_t(BlockSeedsIndex) vec_bseeds_index;

#endif // WORD_FINDER_AUX_H
