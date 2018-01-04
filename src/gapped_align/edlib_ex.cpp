#include "edlib_ex.h"

#include <cassert>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <iostream>

#include "../common/aux_tools.h"
#include "../common/smart_assert.h"

using namespace std;

namespace ns_edlib_ex {

#ifndef GAP_CHAR
#define GAP_CHAR ('-')
#endif

// Status codes
#define EDLIB_STATUS_OK 0
#define EDLIB_STATUS_ERROR 1

// Edit operations.
#define EDLIB_EDOP_MATCH 0    //!< Match.
#define EDLIB_EDOP_INSERT 1   //!< Insertion to target = deletion from query.
#define EDLIB_EDOP_DELETE 2   //!< Deletion from target = insertion to query.
#define EDLIB_EDOP_MISMATCH 3 //!< Mismatch.

inline int
calc_num_blocks(const int n)
{
    return (n + WORD_SIZE - 1) / WORD_SIZE;
}

inline void
calc_block_cell_scores(const Block block, int scores[])
{
    int score = block.score;
    Word mask = HIGH_BIT_MASK;
    for (int i = 0; i < WORD_SIZE - 1; ++i) {
        scores[i] = score;
        if (block.P & mask) --score;
        if (block.M & mask) ++score;
        mask >>= 1;
    }
    scores[WORD_SIZE - 1] = score;
}

void
build_peq(const char* query, const int query_size, Word* peq)
{
    const int nblk = calc_num_blocks(query_size);
    for (int s = 0; s <= AlphabetSize; ++s) {
        for (int b = 0; b < nblk; ++b) {
            const int bid = s * MaxNumBlocks + b;
            if (s < AlphabetSize) {
                peq[bid] = 0;
                for (int r = (b + 1) * WORD_SIZE - 1; r >= b * WORD_SIZE; --r) {
                    peq[bid] <<= 1;
                    if (r >= query_size || query[r] == s) peq[bid] += 1;
                }
            } else {
                peq[bid] = (Word)-1;
            }
        }
    }
}

/**
 * Corresponds to Advance_Block function from Myers.
 * Calculates one word(block), which is part of a column.
 * Highest bit of word (one most to the left) is most bottom cell of block from column.
 * Pv[i] and Mv[i] define vin of cell[i]: vin = cell[i] - cell[i-1].
 * @param [in] Pv  Bitset, Pv[i] == 1 if vin is +1, otherwise Pv[i] == 0.
 * @param [in] Mv  Bitset, Mv[i] == 1 if vin is -1, otherwise Mv[i] == 0.
 * @param [in] Eq  Bitset, Eq[i] == 1 if match, 0 if mismatch.
 * @param [in] hin  Will be +1, 0 or -1.
 * @param [out] PvOut  Bitset, PvOut[i] == 1 if vout is +1, otherwise PvOut[i] == 0.
 * @param [out] MvOut  Bitset, MvOut[i] == 1 if vout is -1, otherwise MvOut[i] == 0.
 * @param [out] hout  Will be +1, 0 or -1.
 */
static inline int calculateBlock(Word Pv, Word Mv, Word Eq, const int hin,
                                 Word &PvOut, Word &MvOut) {
    // hin can be 1, -1 or 0.
    // 1  -> 00...01
    // 0  -> 00...00
    // -1 -> 11...11 (2-complement)

    Word hinIsNeg = (Word)(hin >> 2) & WORD_1; // 00...001 if hin is -1, 00...000 if 0 or 1

    Word Xv = Eq | Mv;
    // This is instruction below written using 'if': if (hin < 0) Eq |= (Word)1;
    Eq |= hinIsNeg;
    Word Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;

    Word Ph = Mv | ~(Xh | Pv);
    Word Mh = Pv & Xh;

    int hout = 0;
    // This is instruction below written using 'if': if (Ph & HIGH_BIT_MASK) hout = 1;
    hout = (Ph & HIGH_BIT_MASK) >> (WORD_SIZE - 1);
    // This is instruction below written using 'if': if (Mh & HIGH_BIT_MASK) hout = -1;
    hout -= (Mh & HIGH_BIT_MASK) >> (WORD_SIZE - 1);

    Ph <<= 1;
    Mh <<= 1;

    // This is instruction below written using 'if': if (hin < 0) Mh |= (Word)1;
    Mh |= hinIsNeg;
    // This is instruction below written using 'if': if (hin > 0) Ph |= (Word)1;
    Ph |= (Word)((hin + 1) >> 1);

    PvOut = Mh | ~(Xv | Ph);
    MvOut = Ph & Xv;

    return hout;
}

void
calc_edit_distance_semi_global(const char* query, const int query_size,
							   const char* target, const int target_size,
							   const Word* peq,
							   int k,
							   const EdlibAlignMode mode,
							   Block* blocks,
							   int& edit_distance,
							   vector<int>& end_positions)
{
    edit_distance = -1;
    end_positions.clear();
    const int nblk = calc_num_blocks(query_size);
    const int W = nblk * WORD_SIZE - query_size;
    int fblk = 0;
    int lblk = min(calc_num_blocks(k + 1), nblk) - 1;
    Block* blk = blocks;

    if (mode == EDLIB_MODE_HW) {
        k = min(query_size, k);
    }

    for (int b = 0; b <= lblk; ++b) {
        blk->score = (b + 1) * WORD_SIZE;
        blk->P = (Word)-1;
        blk->M = (Word)0;
        ++blk;
    }

    int best_d = -1;
    const int start_hout = mode == EDLIB_MODE_HW ? 0 : 1;
    for (int c = 0; c < target_size; ++c) {
        const int tc = target[c];
        const Word* cpeq = peq + tc * MaxNumBlocks;
        int hout = start_hout;
        blk = blocks + fblk;
        cpeq += fblk;
        for (int b = fblk; b <= lblk; ++b) {
            hout = calculateBlock(blk->P, blk->M, *cpeq, hout, blk->P, blk->M);
            blk->score += hout;
            ++blk;
            ++cpeq;
        }
        --blk;
        --cpeq;

        if ( (lblk < nblk - 1) && (blk->score - hout <= k) && ((*(cpeq + 1) & WORD_1) || hout < 0) ) {
            ++lblk;
            ++blk;
            ++cpeq;
            blk->P = (Word)-1;
            blk->M = (Word)0;
            blk->score = (blk - 1)->score - hout + WORD_SIZE + calculateBlock(blk->P, blk->M, *cpeq, hout, blk->P, blk->M);
        } else {
            while (lblk >= fblk && blk->score >= k + WORD_SIZE) {
                --lblk;
                --blk;
                --cpeq;
            }
        }

        if (mode != EDLIB_MODE_HW) {
            while (fblk <= lblk && blocks[fblk].score >= k + WORD_SIZE) {
                ++fblk;
            }
        }

        if (mode == EDLIB_MODE_HW) {
            lblk = max(0, lblk);
        }

        if (lblk < fblk) {
            edit_distance = best_d;
            return;
        }

        if (lblk == nblk - 1) {
            int col_score = blk->score;
            if (col_score <= k) {
                if (best_d == -1 || col_score <= best_d) {
                    if (col_score != best_d) {
                        end_positions.clear();
                        best_d = col_score;
                        k = best_d;
                    }
                    end_positions.push_back(c - W);
                }
            }
        }
    }

    if (lblk == nblk - 1) {
        int scores[WORD_SIZE];
        calc_block_cell_scores(*blk, scores);
        for (int i = 0; i < W; ++i) {
            int col_score = scores[i + 1];
            if (col_score <= k && (best_d == -1 || col_score <= best_d)) {
                if (col_score != best_d) {
                    end_positions.clear();
                    k = best_d = col_score;
                }
                end_positions.push_back(target_size - W + i);
            }
        }
    }

    edit_distance = best_d;
}

void
calc_edit_distance_nw(const char* query, const int query_size,
					  const char* target, const int target_size,
					  const Word* peq,
					  int k,
					  const bool traceback,
					  AlignmentData* align_data,
					  Block* blocks,
					  const int target_stop_position,
					  int& edit_distance,
					  int& end_position)
{
    edit_distance = -1;
    end_position = -1;

    if (target_stop_position > -1 && traceback) {
        cerr << "invalid parameters: "
             << "target_stop_position = " << target_stop_position
             << ", traceback = " << traceback
             << endl;
        exit (1);
    }

    if (k < abs(target_size - query_size)) {
        return;
    }

    k = min(k, max(query_size, target_size));

    const int nblk = calc_num_blocks(query_size);
    const int W = nblk * WORD_SIZE - query_size;
    int fblk = 0;
    int lblk = min(nblk, calc_num_blocks(min(k, (k + query_size - target_size) / 2) + 1)) - 1;
    Block* blk = blocks;

    for (int b = 0; b <= lblk; ++b) {
        blk->score = (b + 1) * WORD_SIZE;
        blk->P = (Word)-1;
        blk->M = (Word)0;
        ++blk;
    }

    for (int c = 0; c < target_size; ++c) {
        const int tc = target[c];
        const Word* cpeq = peq + tc * MaxNumBlocks;
        int hout = 1;
        blk = blocks + fblk;
        for (int b = fblk; b <= lblk; ++b) {
            hout = calculateBlock(blk->P, blk->M, cpeq[b], hout, blk->P, blk->M);
            blk->score += hout;
            ++blk;
        }
        --blk;

        k = min(k, blk->score
                + max(target_size -c - 1, query_size - ((1 + lblk) * WORD_SIZE - 1) - 1)
                + (lblk == nblk - 1 ? W : 0));

        if (lblk + 1 < nblk) {
            bool r = (lblk + 1) * WORD_SIZE - 1 > k - blk->score + 2 * WORD_SIZE - 2 - target_size + c + query_size;
            if (!r) {
                ++lblk;
                ++blk;
                blk->P = (Word)-1;
                blk->M = (Word)0;
                int newHout = calculateBlock(blk->P, blk->M, cpeq[lblk], hout, blk->P, blk->M);
                blk->score = (blk - 1)->score - hout + WORD_SIZE + newHout;
                hout = newHout;
            }
        }

        while (lblk >= fblk
				&&
				(blk->score >= k + WORD_SIZE
				 || ((lblk + 1) * WORD_SIZE -1 > 
					 k - blk->score + 2 * WORD_SIZE - 2 - target_size + c + query_size + 1))) {
            --lblk;
            --blk;
        }

        while (fblk <= lblk
				&&
				(blocks[fblk].score >= k + WORD_SIZE
				 || ((fblk + 1) * WORD_SIZE - 1 <
					 blocks[fblk].score - k - target_size + query_size + c))) {
            ++fblk;
        }

        if (lblk < fblk) {
            edit_distance = end_position = -1;
            return;
        }

        if (traceback) {
            blk = blocks + fblk;
            for (int b = fblk; b <= lblk; ++b) {
                align_data->Ps[MaxNumBlocks * c + b] = blk->P;
                align_data->Ms[MaxNumBlocks * c + b] = blk->M;
                align_data->scores[MaxNumBlocks * c + b] = blk->score;
                align_data->first_blocks[c] = fblk;
                align_data->last_blocks[c] = lblk;
                ++blk;
            }
        }

        if (c == target_stop_position) {
            for (int b = fblk; b <= lblk; ++b) {
                align_data->Ps[b] = blocks[b].P;
                align_data->Ms[b] = blocks[b].M;
                align_data->scores[b] = blocks[b].score;
                align_data->first_blocks[0] = fblk;
                align_data->last_blocks[0] = lblk;
            }
            edit_distance = -1;
            end_position = target_stop_position;
            return;
        }
    }

    if (lblk == nblk - 1) {
        int scores[WORD_SIZE];
        calc_block_cell_scores(blocks[lblk], scores);
        int col_score = scores[W];
        if (col_score <= k) {
            edit_distance = col_score;
            end_position = target_size - 1;
            return;
        }
    }

    edit_distance = end_position = -1;
}

/**
 * Finds one possible alignment that gives optimal score by moving back through the dynamic programming matrix,
 * that is stored in alignData. Consumes large amount of memory: O(queryLength * targetLength).
 * @param [in] queryLength  Normal length, without W.
 * @param [in] targetLength  Normal length, without W.
 * @param [in] bestScore  Best score.
 * @param [in] alignData  Data obtained during finding best score that is useful for finding alignment.
 * @param [out] alignment  Alignment.
 * @param [out] alignmentLength  Length of alignment.
 * @return Status code.
 */
static int obtainAlignmentTraceback(const int queryLength, const int targetLength,
                                    const int bestScore, const AlignmentData* const alignData,
                                    vector<unsigned char>& alignment) {
	const int nblk = calc_num_blocks(queryLength);
	const int W = nblk * WORD_SIZE - queryLength;
	alignment.clear();
    int c = targetLength - 1; // index of column
    int b = nblk - 1; // index of block in column
    int currScore = bestScore; // Score of current cell
    int lScore  = -1; // Score of left cell
    int uScore  = -1; // Score of upper cell
    int ulScore = -1; // Score of upper left cell
    Word currP = alignData->Ps[c * MaxNumBlocks + b]; // P of current block
    Word currM = alignData->Ms[c * MaxNumBlocks + b]; // M of current block
    // True if block to left exists and is in band
    bool thereIsLeftBlock = c > 0 && b >= alignData->first_blocks[c-1] && b <= alignData->last_blocks[c-1];
    // We set initial values of lP and lM to 0 only to avoid compiler warnings, they should not affect the
    // calculation as both lP and lM should be initialized at some moment later (but compiler can not
    // detect it since this initialization is guaranteed by "business" logic).
    Word lP = 0, lM = 0;
    if (thereIsLeftBlock) {
        lP = alignData->Ps[(c - 1) * MaxNumBlocks + b]; // P of block to the left
        lM = alignData->Ms[(c - 1) * MaxNumBlocks + b]; // M of block to the left
    }
    currP <<= W;
    currM <<= W;
    int blockPos = WORD_SIZE - W - 1; // 0 based index of current cell in blockPos

    // TODO(martin): refactor this whole piece of code. There are too many if-else statements,
    // it is too easy for a bug to hide and to hard to effectively cover all the edge-cases.
    // We need better separation of logic and responsibilities.
    while (true) {
        if (c == 0) {
            thereIsLeftBlock = true;
            lScore = b * WORD_SIZE + blockPos + 1;
            ulScore = lScore - 1;
        }

        // TODO: improvement: calculate only those cells that are needed,
        //       for example if I calculate upper cell and can move up,
        //       there is no need to calculate left and upper left cell
        //---------- Calculate scores ---------//
        if (lScore == -1 && thereIsLeftBlock) {
            lScore = alignData->scores[(c - 1) * MaxNumBlocks + b]; // score of block to the left
            for (int i = 0; i < WORD_SIZE - blockPos - 1; i++) {
                if (lP & HIGH_BIT_MASK) lScore--;
                if (lM & HIGH_BIT_MASK) lScore++;
                lP <<= 1;
                lM <<= 1;
            }
        }
        if (ulScore == -1) {
            if (lScore != -1) {
                ulScore = lScore;
                if (lP & HIGH_BIT_MASK) ulScore--;
                if (lM & HIGH_BIT_MASK) ulScore++;
            }
            else if (c > 0 && b-1 >= alignData->first_blocks[c-1] && b-1 <= alignData->last_blocks[c-1]) {
                // This is the case when upper left cell is last cell in block,
                // and block to left is not in band so lScore is -1.
                ulScore = alignData->scores[(c - 1) * MaxNumBlocks + b - 1];
            }
        }
        if (uScore == -1) {
            uScore = currScore;
            if (currP & HIGH_BIT_MASK) uScore--;
            if (currM & HIGH_BIT_MASK) uScore++;
            currP <<= 1;
            currM <<= 1;
        }
        //-------------------------------------//

        // TODO: should I check if there is upper block?

        //-------------- Move --------------//
        // Move up - insertion to target - deletion from query
        if (uScore != -1 && uScore + 1 == currScore) {
            currScore = uScore;
            lScore = ulScore;
            uScore = ulScore = -1;
            if (blockPos == 0) { // If entering new (upper) block
                if (b == 0) { // If there are no cells above (only boundary cells)
					alignment.push_back(EDLIB_EDOP_INSERT);
					for (int i = 0; i < c + 1; ++i) alignment.push_back(EDLIB_EDOP_DELETE);
					break;
                } else {
                    blockPos = WORD_SIZE - 1;
                    b--;
                    currP = alignData->Ps[c * MaxNumBlocks + b];
                    currM = alignData->Ms[c * MaxNumBlocks + b];
                    if (c > 0 && b >= alignData->first_blocks[c-1] && b <= alignData->last_blocks[c-1]) {
                        thereIsLeftBlock = true;
                        lP = alignData->Ps[(c - 1) * MaxNumBlocks + b]; // TODO: improve this, too many operations
                        lM = alignData->Ms[(c - 1) * MaxNumBlocks + b];
                    } else {
                        thereIsLeftBlock = false;
                        // TODO(martin): There may not be left block, but there can be left boundary - do we
                        // handle this correctly then? Are l and ul score set correctly? I should check that / refactor this.
                    }
                }
            } else {
                blockPos--;
                lP <<= 1;
                lM <<= 1;
            }
            // Mark move
            ///-(*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
			alignment.push_back(EDLIB_EDOP_INSERT);
        }
        // Move left - deletion from target - insertion to query
        else if (lScore != -1 && lScore + 1 == currScore) {
            currScore = lScore;
            uScore = ulScore;
            lScore = ulScore = -1;
            c--;
            if (c == -1) { // If there are no cells to the left (only boundary cells)
                ///(*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE; // Move left
                ///int numUp = b * WORD_SIZE + blockPos + 1;
                ///for (int i = 0; i < numUp; i++) // Move up until end
                ///    (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
                ///break;
				
				alignment.push_back(EDLIB_EDOP_DELETE);
				int numUp = b * WORD_SIZE + blockPos + 1;
				for (int i = 0; i < numUp; i++) alignment.push_back(EDLIB_EDOP_INSERT);
				break;
            }
            currP = lP;
            currM = lM;
            if (c > 0 && b >= alignData->first_blocks[c-1] && b <= alignData->last_blocks[c-1]) {
                thereIsLeftBlock = true;
                lP = alignData->Ps[(c - 1) * MaxNumBlocks + b];
                lM = alignData->Ms[(c - 1) * MaxNumBlocks + b];
            } else {
                if (c == 0) { // If there are no cells to the left (only boundary cells)
                    thereIsLeftBlock = true;
                    lScore = b * WORD_SIZE + blockPos + 1;
                    ulScore = lScore - 1;
                } else {
                    thereIsLeftBlock = false;
                }
            }
            // Mark move
            ///(*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE;
			alignment.push_back(EDLIB_EDOP_DELETE);
        }
        // Move up left - (mis)match
        else if (ulScore != -1) {
            unsigned char moveCode = ulScore == currScore ? EDLIB_EDOP_MATCH : EDLIB_EDOP_MISMATCH;
            currScore = ulScore;
            uScore = lScore = ulScore = -1;
            c--;
            if (c == -1) { // If there are no cells to the left (only boundary cells)
                ///(*alignment)[(*alignmentLength)++] = moveCode; // Move left
                ///int numUp = b * WORD_SIZE + blockPos;
                ///for (int i = 0; i < numUp; i++) // Move up until end
                ///    (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
                ///break;
				
				alignment.push_back(moveCode);
				int numUp = b * WORD_SIZE + blockPos;
				for (int i = 0; i < numUp; i++) alignment.push_back(EDLIB_EDOP_INSERT);
				break;
            }
            if (blockPos == 0) { // If entering upper left block
                if (b == 0) { // If there are no more cells above (only boundary cells)
                    ///(*alignment)[(*alignmentLength)++] = moveCode; // Move up left
                    ///for (int i = 0; i < c + 1; i++) // Move left until end
                    ///    (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE;
                    ///break;
					
					alignment.push_back(moveCode);
					for (int i = 0; i < c + 1; i++) alignment.push_back(EDLIB_EDOP_DELETE);
					break;
                }
                blockPos = WORD_SIZE - 1;
                b--;
                currP = alignData->Ps[c * MaxNumBlocks + b];
                currM = alignData->Ms[c * MaxNumBlocks + b];
            } else { // If entering left block
                blockPos--;
                currP = lP;
                currM = lM;
                currP <<= 1;
                currM <<= 1;
            }
            // Set new left block
            if (c > 0 && b >= alignData->first_blocks[c-1] && b <= alignData->last_blocks[c-1]) {
                thereIsLeftBlock = true;
                lP = alignData->Ps[(c - 1) * MaxNumBlocks + b];
                lM = alignData->Ms[(c - 1) * MaxNumBlocks + b];
            } else {
                if (c == 0) { // If there are no cells to the left (only boundary cells)
                    thereIsLeftBlock = true;
                    lScore = b * WORD_SIZE + blockPos + 1;
                    ulScore = lScore - 1;
                } else {
                    thereIsLeftBlock = false;
                }
            }
            // Mark move
            ///(*alignment)[(*alignmentLength)++] = moveCode;
			alignment.push_back(moveCode);
        } else {
            // Reached end - finished!
            break;
        }
        //----------------------------------//
    }

    reverse(alignment.begin(), alignment.end());
    return EDLIB_STATUS_OK;
}

void
edlibAlignmentToCigar1(const unsigned char* alignment, 
					   int alignmentLength, 
					   EdlibCigarFormat cigarFormat,
					   vector<GapAlignOp>& cigar)
{
	cigar.clear();
	if (cigarFormat != EDLIB_CIGAR_EXTENDED && cigarFormat != EDLIB_CIGAR_STANDARD) {
        return;
    }

    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = {'=', 'I', 'D', 'X'};
    if (cigarFormat == EDLIB_CIGAR_STANDARD) {
        moveCodeToChar[0] = moveCodeToChar[3] = 'M';
    }

    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
	GapAlignOp op;
    for (int i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if (i == alignmentLength || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0)) {
			op.num = numOfSameMoves;
			op.op = lastMove;
			cigar.push_back(op);
            // If not at the end, start new sequence of moves.
            if (i < alignmentLength) {
                // Check if alignment has valid values.
                if (alignment[i] > 3) {
                    cigar.clear();
                    return;
                }
                numOfSameMoves = 0;
            }
        }
        if (i < alignmentLength) {
            lastMove = moveCodeToChar[alignment[i]];
            numOfSameMoves++;
        }
    }
}

bool 
EdLibAligner::align(const char* query, const int query_size,
				   const char* target, const int target_size,
				   const double error,
				   char* query_align,
				   char* target_align,
				   int& qend,
				   int& tend)
{
	return go(query, query_size, target, target_size, error, query_align, target_align, qend, tend);
}

bool 
EdLibAligner::go(const char* query, const int query_size,
				const char* target, const int target_size,
				const double error,
				char* query_align,
				char* target_align,
				int& qend,
				int& tend)
{
	query_align[0] = '\0';
	target_align[0] = '\0';
	qend = 0;
	tend = 0;
	
	result.clear();
	build_peq(query, query_size, peq);
	int k = min(query_size, target_size) * error * 1.1;
	calc_edit_distance_semi_global(query, query_size,
								   target, target_size,
								   peq,
								   k,
								   EDLIB_MODE_SHW,
								   blocks,
								   result.edit_distance,
								   result.end_locations);
	if (result.edit_distance == -1) return 0;
	assert(!result.end_locations.empty());
	
	int edit_distance;
	int end_position;
	calc_edit_distance_nw(query, query_size,
						  target, result.end_locations[0] + 1,
						  peq,
						  result.edit_distance,
						  1,
						  &align_data,
						  blocks,
						  -1,
						  edit_distance,
						  end_position);
	assert(edit_distance == result.edit_distance);
	assert(end_position == result.end_locations[0]);
	
	obtainAlignmentTraceback(query_size, 
							 result.end_locations[0] + 1,
							 edit_distance, 
							 &align_data,
							 result.alignment);
	edlibAlignmentToCigar1(result.alignment.data(), result.alignment.size(), EDLIB_CIGAR_STANDARD, cigar);
	const char* q = query;
	const char* t = target;
	int aidx = 0;
	for (size_t i = 0; i != cigar.size(); ++i) {
		int n = cigar[i].num;
		char op = cigar[i].op;
		switch (op) {
			case 'M':
				for (int c = 0; c < n; ++c) {
					int qc = *q++;
					int tc = *t++;
					assert(qc >= 0 && qc < 4);
					assert(tc >= 0 && tc < 4);
					query_align[aidx] = DNADecoder::decode(qc); //dt[qc];
					target_align[aidx] = DNADecoder::decode(tc); //dt[tc];
					++qend;
					++tend;
					++aidx;
				}
				break;
			case 'I':
				for (int c = 0; c < n; ++c) {
					int qc = *q++;
					assert(qc >= 0 && qc < 4);
					query_align[aidx] = DNADecoder::decode(qc); //dt[qc];
					target_align[aidx] = GAP_CHAR;
					++aidx;
					++qend;
				}
				break;
			case 'D':
				for (int c = 0; c < n; ++c) {
					int tc = *t++;
					assert(tc >= 0 && tc < 4);
					query_align[aidx] = GAP_CHAR;
					target_align[aidx] = DNADecoder::decode(tc); //dt[tc];
					++aidx;
					++tend;
				}
				break;
			default:
				cout << n << '\t' << (int)op << "\n";
				assert(0);
				break;
		}
	}
	query_align[aidx] = '\0';
	target_align[aidx] = '\0';
	
	return 1;
}

} // ns_edlib_ex
