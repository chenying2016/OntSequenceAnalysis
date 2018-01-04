#ifndef EDLIB_EX_H
#define EDLIB_EX_H

#include <stdint.h>
#include <string>
#include <vector>

namespace ns_edlib_ex {

typedef uint64_t Word;
#define WORD_SIZE 64
#define WORD_1 ((Word)1)
#define HIGH_BIT_MASK (WORD_1 << (WORD_SIZE - 1))

///*
#define MaxSeqSize 1024
#define MaxNumBlocks 16
#define AlphabetSize 4
//*/

/*
#define MaxSeqSize 2048
#define MaxNumBlocks 32
#define AlphabetSize 4
*/

/**
 * Describes cigar format.
 * @see http://samtools.github.io/hts-specs/SAMv1.pdf
 * @see http://drive5.com/usearch/manual/cigar.html
 */
enum EdlibCigarFormat{
	EDLIB_CIGAR_STANDARD,  //!< Match: 'M', Insertion: 'I', Deletion: 'D', Mismatch: 'M'.
	EDLIB_CIGAR_EXTENDED   //!< Match: '=', Insertion: 'I', Deletion: 'D', Mismatch: 'X'.
};

enum EdlibAlignMode
{
    EDLIB_MODE_NW,
    EDLIB_MODE_SHW,
    EDLIB_MODE_HW
};

struct AlignmentData
{
    Word* Ps;
    Word* Ms;
    int* scores;
    int* first_blocks;
    int* last_blocks;

    AlignmentData() {
        Ps = new Word[MaxNumBlocks * MaxSeqSize];
        Ms = new Word[MaxNumBlocks * MaxSeqSize];
        scores = new int[MaxNumBlocks * MaxSeqSize];
        first_blocks = new int[MaxSeqSize];
        last_blocks = new int[MaxSeqSize];
    }

    ~AlignmentData() {
        delete[] Ps;
        delete[] Ms;
        delete[] scores;
        delete[] first_blocks;
        delete[] last_blocks;
    }
};

struct Block
{
    Word P;
    Word M;
    int score;
};

struct AlignResult
{
    int edit_distance;
    std::vector<int> end_locations;
    std::vector<int> start_locations;
    std::vector<unsigned char> alignment;

    void clear() {
        edit_distance = -1;
        end_locations.clear();
        start_locations.clear();
        alignment.clear();
    }
};

#define OP_SUB ('M')
#define OP_INS ('I')
#define OP_DEL ('D')
struct GapAlignOp
{
	int num;
	int op;
};

class EdLibAligner
{
public:
	EdLibAligner() {
		blocks = new Block[MaxNumBlocks];
		peq = new Word[(AlphabetSize + 1) * MaxNumBlocks];
	}
	~EdLibAligner() {
		delete[] blocks;
		delete[] peq;
	}
	
	bool align(const char* query, const int query_size,
			   const char* target, const int target_size,
			   const double error,
			   char* query_align,
			   char* target_align,
			   int& qend,
			   int& tend);
	
private:
	bool go(const char* query, const int query_size,
			const char* target, const int target_size,
			const double error,
			char* query_align,
			char* target_align,
		    int& qend,
		    int& tend);
	
private:
	AlignmentData 	align_data;
	AlignResult		result;
	Block*			blocks;
	Word*			peq;
	std::vector<GapAlignOp> cigar;
};

} // ns_edlib_ex

using ns_edlib_ex::EdLibAligner;

#endif // EDLIB_EX_H
