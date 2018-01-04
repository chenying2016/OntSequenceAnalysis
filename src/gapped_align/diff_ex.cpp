#include "diff_ex.h"
#include "../common/smart_assert.h"

using namespace std;

namespace ns_diff_ex
BEGIN_NAMESPACE_SCOPE

struct SCompareDPathData2
{
    bool operator()(const DPathData2& a, const DPathData2& b) {
        return (a.d == b.d) ? (a.k < b.k) : (a.d < b.d);
    }
};

struct DPathDataExtractor
{
	DPathData2* operator()(const int d, const int k)
	{
		target.d = d;
		target.k = k;
		return std::lower_bound(d_path_list, d_path_list + d_path_list_size, target, scmp);
	}
	
	DPathDataExtractor(DPathData2* d_path, const int n)
		: d_path_list(d_path),
		  d_path_list_size(n) {
	}
	
    DPathData2 target;
	SCompareDPathData2 scmp;
	DPathData2* d_path_list;
	const int d_path_list_size;
};

template <typename T>
inline T extract_char(const char* A, int i, bool forward)
{
    if (forward) {
        return static_cast<T>(A[i]);
    } else {
        return static_cast<T>(A[-i]);
    }
}

void
GetAlignString(const char* query, const int q_len, 
			   const char* target, const int t_len,
			   DPathData2* d_path, const int d_path_size,
			   PathPoint* aln_path, Alignment* align,
			   const int qe, const int te, 
			   int d, int k,
			   const int right_extend)
{
	int cd = d;
	int ck = k;
	int aln_path_idx = 0;
	int i;
	DPathDataExtractor dpath_extrator(d_path, d_path_size);
	while (cd >= 0 && aln_path_idx < q_len + t_len + 1)
	{
		DPathData2* d_path_aux = dpath_extrator(cd, ck);
		aln_path[aln_path_idx].x = d_path_aux->x2;
		aln_path[aln_path_idx].y = d_path_aux->y2;
		++aln_path_idx;
		aln_path[aln_path_idx].x = d_path_aux->x1;
		aln_path[aln_path_idx].y = d_path_aux->y1;
		++aln_path_idx;
		ck = d_path_aux->pre_k;
		cd -= 1;
	}
	--aln_path_idx;
	int cx = aln_path[aln_path_idx].x;
	int cy = aln_path[aln_path_idx].y;
	align->aln_q_s = cx;
	align->aln_t_s = cy;
	int aln_pos = 0;
	while (aln_path_idx > 0)
	{
		--aln_path_idx;
		int nx = aln_path[aln_path_idx].x;
		int ny = aln_path[aln_path_idx].y;
		if (cx == nx && cy == ny) continue;
		
		if (cx == nx && cy != ny) {
			for (i = 0; i < ny - cy; ++i) {
				align->q_aln_str[aln_pos + i] = GAP_CODE;
				align->t_aln_str[aln_pos + i] = extract_char<char>(target, cy + i, right_extend);
			}
			aln_pos += ny - cy;
		} else if (cx != nx && cy == ny) {
			for (i = 0; i < nx - cx; ++i) {
				align->q_aln_str[aln_pos + i] = extract_char<char>(query, cx + i, right_extend);
				align->t_aln_str[aln_pos + i] = GAP_CODE;
			}
			aln_pos += nx - cx;
		} else {
			for (i = 0; i < nx - cx; ++i) {
				align->q_aln_str[aln_pos + i] = extract_char<char>(query, cx + i, right_extend);
			}
			for (i = 0; i < ny - cy; ++i) {
				align->t_aln_str[aln_pos + i] = extract_char<char>(target, cy + i, right_extend);
			}
			aln_pos += ny - cy;
		}
		
		cx = nx;
		cy = ny;
	}
	align->aln_str_size = aln_pos;
}


int Align(const char* query, const int q_len, const char* target, const int t_len,
		  const int band_tolerance, const int get_aln_str, Alignment* align,
		  int* V, int* U, DPathData2* d_path, PathPoint* aln_path, const int right_extend)
{
	int k_offset;
	int  d;
	int  k, k2;
	int best_m;
	int min_k, new_min_k, max_k, new_max_k, pre_k;
	int x = -1, y = -1;
	int max_d, band_size;
	unsigned long d_path_idx = 0, max_idx = 0;
	int aligned = 0;
	int best_x = -1, best_y = -1, best_d = q_len + t_len + 100, best_k = 0, best_d_path_idx = -1;

	max_d = (int)(.4 * (q_len + t_len));
	k_offset = max_d;
	band_size = band_tolerance * 2;
	align->init();
	best_m = -1;
	min_k = 0;
	max_k = 0;
	d_path_idx = 0;
	max_idx = 0;

	for (d = 0; d < max_d; ++d)
	{
		if (max_k - min_k > band_size) break;
		
		for (k = min_k; k <= max_k; k += 2)
		{
			if( k == min_k || (k != max_k && V[k - 1 + k_offset] < V[k + 1 + k_offset]) )
			{ pre_k = k + 1; x = V[k + 1 + k_offset]; }
			else
			{ pre_k = k - 1; x = V[k - 1 + k_offset] + 1; }
			y = x - k;
			d_path[d_path_idx].d = d;
			d_path[d_path_idx].k = k;
			d_path[d_path_idx].x1 = x;
			d_path[d_path_idx].y1 = y;

			if (right_extend) {
				while( x < q_len && y < t_len && query[x] == target[y]) { ++x; ++y; }
			} else {
				while( x < q_len && y < t_len && query[-x] == target[-y]) { ++x; ++y; }
			}

			d_path[d_path_idx].x2 = x;
			d_path[d_path_idx].y2 = y;
			d_path[d_path_idx].pre_k = pre_k;
			++d_path_idx;

			V[k + k_offset] = x;
			U[k + k_offset] = x + y;
			if (x + y > best_m) {
				best_m = x + y;
				best_x = x;
				best_y = y;
				best_d = d;
				best_k = k;
				best_d_path_idx = d_path_idx;
			}
			if (x >= q_len || y >= t_len)
			{ aligned = 1; max_idx = d_path_idx; break; }
		}

		// for banding
		new_min_k = max_k;
		new_max_k = min_k;
		for (k2 = min_k; k2 <= max_k; k2 += 2)
			if (U[k2 + k_offset] >= best_m - band_tolerance)
			{ new_min_k = std::min(new_min_k, k2); new_max_k = std::max(new_max_k, k2); }
		max_k = new_max_k + 1;
		min_k = new_min_k - 1;

		if (aligned)
		{
			align->aln_q_e = x;
			align->aln_t_e = y;
			align->dist = d;
			align->aln_str_size = (x + y + d) / 2;
			align->aln_q_s = 0;
			align->aln_t_s = 0;

			if (get_aln_str) {
				GetAlignString(query, q_len, target, t_len, d_path, max_idx, aln_path, align, x, y, d, k, right_extend);
			}
			break;
		} 
	}
	
	if (!aligned) {
		if (best_x > 0) {
			align->aln_q_e = best_x;
			align->aln_t_e = best_y;
			align->dist = best_d;
			align->aln_str_size = (best_x + best_y + best_d) / 2;
			align->aln_q_s = 0;
			align->aln_t_s = 0;
			if (get_aln_str) {
				GetAlignString(query, q_len, target, t_len, d_path, best_d_path_idx, aln_path, align, best_x, best_y, best_d, best_k, right_extend);
			}
		} else {
			align->aln_q_e = 0;
			align->aln_t_e = 0;
			align->dist = 0;
			align->aln_str_size = 0;
			align->aln_q_s = 0;
			align->aln_t_s = 0;
		}
	}
	
	return (align->aln_q_e == q_len || align->aln_t_e == t_len);
}

bool 
DiffAligner::align(const char* query, 
		   const int query_size,
		   const char* target, 
		   const int target_size,
		   const double error,
		   char* query_align,
		   char* target_align,
		   int& qend,
		   int& tend)
{
	clean();
	int band_tolerance = max(0.3, error) * max(target_size, query_size);
	Align(query, 
		  query_size, 
		  target, 
		  target_size,
		  band_tolerance,
		  1,
		  result,
		  dynq,
		  dynt,
		  d_path,
		  aln_path,
		  1);
	qend = result->aln_q_e;
	tend = result->aln_t_e;
	int i;
	for (i = 0; i < result->aln_str_size; ++i) {
		query_align[i] = DNADecoder::decode(result->q_aln_str[i]);
		target_align[i] = DNADecoder::decode(result->t_aln_str[i]);
	}
	query_align[i] = '\0';
	target_align[i] = '\0';
	
	{
		int qid = 0, tid = 0;
		for (i = 0; i < result->aln_str_size; ++i) {
			int qc = result->q_aln_str[i];
			int tc = result->t_aln_str[i];
			if (qc != GAP_CODE) {
				int qc1 = query[qid];
				r_assert(qc == qc1)(qid)(i)(qc)(qc1)(qend)(tend)(result->aln_str_size);
				++qid;
			}
			if (tc != GAP_CODE) {
				int tc1 = target[tid];
				r_assert(tc == tc1)(tid)(i)(tc)(tc1)(qend)(tend)(result->aln_str_size);
				++tid;
			}
		}
		r_assert(qid == qend)(qid)(qend)(query_size)(target_size);
		r_assert(tid == tend)(tid)(tend)(query_size)(target_size);
	}
	
	return 1;
}

END_NAMESPACE_SCOPE
