#ifndef DIFF_EX_H
#define DIFF_EX_H

#include "../common/aux_tools.h"

#include <string>

namespace ns_diff_ex
BEGIN_NAMESPACE_SCOPE

struct DPathData
{
    int pre_k, x1, y1, x2, y2;
};

struct DPathData2
{
	int d, k, pre_k, x1, y1, x2, y2;
};

struct PathPoint
{
	int x, y;
};

struct Alignment
{
    int aln_str_size;
    int dist;
    int aln_q_s;
    int aln_q_e;
    int aln_t_s;
    int aln_t_e;
    char* q_aln_str;
    char* t_aln_str;

    void init() {
        aln_str_size = 0;
        aln_q_s = aln_q_e = 0;
        aln_t_s = aln_t_e = 0;
    }
	
	Alignment(const idx max_aln_size) {
		snew(q_aln_str, char, max_aln_size);
		snew(t_aln_str, char, max_aln_size);
	}
	
	~Alignment() {
		sfree(q_aln_str);
		sfree(t_aln_str);
	}
};

class DiffAligner
{
public:
	DiffAligner() {
		result = new Alignment(kNumDpRows);
		snew(dynq, int, kNumDpRows);
		snew(dynt, int, kNumDpRows);
		snew(d_path, DPathData2, kNumDpCells);
		snew(aln_path, PathPoint, kNumDpCells);
	}
	
	~DiffAligner() {
		delete result;
		sfree(dynq);
		sfree(dynt);
		sfree(d_path);
		sfree(aln_path);
	}
	
	bool align(const char* query, 
			   const int query_size,
			   const char* target, 
			   const int target_size,
			   const double error,
			   char* query_align,
			   char* target_align,
			   int& qend,
			   int& tend);
	
	void clean() {
		std::fill(dynq, dynq + kNumDpRows, 0);
		std::fill(dynt, dynt + kNumDpRows, 0);
	}
	
private:
	Alignment*				result;
    int*           		 	dynq;
    int*            		dynt;
    DPathData2*     		d_path;
    PathPoint*      		aln_path;
	
	static const std::size_t kNumDpRows = 4096;
	static const std::size_t kNumDpCells = 5000000;
};

END_NAMESPACE_SCOPE

typedef ns_diff_ex::DiffAligner DiffAligner;

#endif // DIFF_EX_H
