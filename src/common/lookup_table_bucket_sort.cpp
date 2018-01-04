#include "lookup_table_bucket_sort.h"

#include "defs.h"
#include "packed_db.h"
#include "timer.h"

#include <algorithm>

using namespace std;

namespace ns_mc_lookup_table {

/// u64 consists of two parts
/// left most 30 bits represent hash value
/// right most 34 bits represent offset

#define HASH_BITS 				30
#define OFFSET_BITS				34
#define HASH_MASK				((((u64)1) << HASH_BITS) - 1)
#define OFFSET_MASK				((((u64)1) << OFFSET_BITS) - 1)
#define EXTRACT_HASH(u)			(((u) >> OFFSET_BITS) & HASH_MASK)
#define EXTRACT_OFFSET(u)		((u) & OFFSET_MASK)
#define PACK_LIST_ITEM(h, o)	(((h) << OFFSET_BITS) | (o))


/// bucket sort parameter
/// BS = bucket sort
static const u64 BS_B = 8;
static const u64 BS_MASK = 255;
static const u64 BS_PASS = 4; // 32/8 = 4
static const u64 BS_BucketSize = 256;

struct BucketSortData
{
	u64* 				src;
	u64* 				trg;
	u64					src_size;
	u64*  				pass;
	int* 				gtid;
	pthread_mutex_t* 	gtid_lock;
	int					num_threads;
	u64**				buckets;
	u64***				next_buckets;
};

#define get_next_hash { \
	const u8 c = ref->get_char(offset + j); \
	hash = ((hash << 2) | c) & hash_mask; \
}

short*
get_kmer_counts(const PackedDB* ref, const int kmer_size, const int kmer_cnt_cutoff, u64& num_kmers)
{
	DynamicTimer dt(__func__);
	
	short* cnt_table = new short[kmer_cnt_cutoff + 2];
	for (short i = 0; i <= kmer_cnt_cutoff; ++i) cnt_table[i] = i + 1;
	cnt_table[kmer_cnt_cutoff + 1] = kmer_cnt_cutoff + 1;
	u64 max_hash = 1; max_hash = max_hash << (kmer_size * 2);
	u64 hash_mask = max_hash - 1;
	short* kmer_cnts = new short[max_hash];
	fill(kmer_cnts, kmer_cnts + max_hash, 0);
	
	const idx ns = ref->num_seqs();
	for (idx i = 0; i < ns; ++i) {
		const idx offset = ref->seq_offset(i);
		const idx size = ref->seq_size(i);
		u64 hash = 0;
		for (idx j = 0; j < kmer_size - 1; ++j) {
			get_next_hash;
		}
		for (idx j = kmer_size - 1; j < size; ++j) {
			get_next_hash;
			kmer_cnts[hash] = cnt_table[kmer_cnts[hash]];
		}
	}
	
	num_kmers = 0;
	for (u64 i = 0; i != max_hash; ++i) {
		if (kmer_cnts[i] > kmer_cnt_cutoff) kmer_cnts[i] = 0;
		num_kmers += kmer_cnts[i];
	}
	
	delete[] cnt_table;
	return kmer_cnts;
}

u64*
get_offset_list(const PackedDB*ref, short* kmer_cnts, const int kmer_size, const u64 num_kmers)
{
	DynamicTimer dt(__func__);
	
	u64* offset_list = new u64[num_kmers];
	u64 cnt = 0;
	u64 max_hash = 1; max_hash = max_hash << (kmer_size * 2);
	u64 hash_mask = max_hash - 1;
	
	const idx ns = ref->num_seqs();
	for (idx i = 0; i < ns; ++i) {
		const idx offset = ref->seq_offset(i);
		const idx size = ref->seq_size(i);
		u64 hash = 0;
		for (idx j = 0; j < kmer_size - 1; ++j) {
			get_next_hash;
		}
		for (idx j = kmer_size - 1; j < size; ++j) {
			get_next_hash;
			u64 k = offset + j + 1 - kmer_size;
			if (kmer_cnts[hash]) offset_list[cnt++] = PACK_LIST_ITEM(hash, k);
		}
	}
	
	return offset_list;
}

/// bucket sort routines

void*
sort_thread(void* arg)
{
	BucketSortData* bsd = (BucketSortData*)(arg);
	int tid;
	pthread_mutex_lock(bsd->gtid_lock);
	tid = *bsd->gtid;
	++(*bsd->gtid);
	pthread_mutex_unlock(bsd->gtid_lock);
	r_assert(tid < bsd->num_threads);
	const u64 shift = BS_B * (*bsd->pass);
	const u64 part = (bsd->src_size + bsd->num_threads - 1) / bsd->num_threads;
	const u64 from = tid * part;
	const u64 to = min(from + part, bsd->src_size);
	u64* src = bsd->src;
	u64* trg = bsd->trg;
	for (u64 i = from; i < to; ++i) {
		u64 c = EXTRACT_HASH(src[i]);
		u64 b = c >> shift;
		u64 b_idx = b & BS_MASK;
		u64 x = bsd->buckets[tid][b_idx]++;
		r_assert(x < bsd->src_size)(tid)(i)(c)(b)(b_idx)(x)(bsd->src_size);
		trg[x] = src[i];
		++bsd->next_buckets[tid][x/part][(b>>BS_B)&BS_MASK];
	}
	
	return NULL;
}

void
init_buckets(u64** buckets, const int num_threads, u64* src, const u64 src_size)
{
	DynamicTimer dt(__func__);
	
	const u64 part = (src_size + num_threads - 1) / num_threads;
	for (int i = 0; i < num_threads; ++i) {
		u64 from = part * i;
		u64 to = min(from + part, src_size);
		for (u64 k = from; k < to; ++k) {
			u64 h = EXTRACT_HASH(src[k]);
			u64 b = h & BS_MASK;
			++buckets[i][b];
		}
	}
}

struct InitBucketsData
{
	u64* 				src;
	u64 				src_size;
	u64** 				buckets;
	int* 				thread_id;
	pthread_mutex_t* 	tid_lock;
	int					num_threads;
};

void*
init_buckets_worker(void* arg)
{
	InitBucketsData* ibd = (InitBucketsData*)(arg);
	int tid;
	pthread_mutex_lock(ibd->tid_lock);
	tid = *ibd->thread_id;
	++(*ibd->thread_id);
	pthread_mutex_unlock(ibd->tid_lock);
	r_assert(tid < ibd->num_threads);
	
	const u64 part = (ibd->src_size + ibd->num_threads - 1) / ibd->num_threads;
	const u64 from = part * tid;
	const u64 to = min(from + part, ibd->src_size);
	for (u64 k = from; k < to; ++k) {
		u64 h = EXTRACT_HASH(ibd->src[k]);
		u64 b = h & BS_MASK;
		++ibd->buckets[tid][b];
	}
	
	return NULL;
}

void
init_buckets_mt(u64** buckets, const int num_threads, u64* src, const u64 src_size)
{
	DynamicTimer dt(__func__);
	
	InitBucketsData ibds[num_threads];
	int thread_id = 0;
	pthread_mutex_t tid_lock;
	pthread_mutex_init(&tid_lock, NULL);
	pthread_t tids[num_threads];
	
	for (int i = 0; i < num_threads; ++i) {
		ibds[i].src 		= src;
		ibds[i].src_size 	= src_size;
		ibds[i].buckets 	= buckets;
		ibds[i].thread_id 	= &thread_id;
		ibds[i].tid_lock 	= &tid_lock;
		ibds[i].num_threads = num_threads;
		pthread_create(tids + i, NULL, init_buckets_worker, ibds + i);
	}
	for (int i = 0; i < num_threads; ++i) {
		pthread_join(tids[i], NULL);
	}
}

void
fill_bucket_from_next(u64** buckets, u64*** next_buckets, const int num_threads)
{
	for (int i = 0; i < num_threads; ++i) {
		for (u64 j = 0; j < BS_BucketSize; ++j) {
			int c = 0;
			for (int k = 0; k < num_threads; ++k) {
				c += next_buckets[k][i][j];
			}
			buckets[i][j] = c;
		}
	}
}

void
update_bucket(u64** buckets, u64** tmp_buckets, const u64 src_size, const int num_threads)
{
	for (int i = 0; i < num_threads; ++i) {
		copy(buckets[i], buckets[i] + BS_BucketSize, tmp_buckets[i]);
	}
	
	for (int t = 0; t < num_threads; ++t) {
		for (u64 b = 0; b < BS_BucketSize; ++b) {
			int cnt = 0;
			for (int u = 0; u < num_threads; ++u) {
				for (u64 c = 0; c < b; ++c) {
					cnt += tmp_buckets[u][c];
				}
			}
			for (int u = 0; u < t; ++u) {
				cnt += tmp_buckets[u][b];
			}
			r_assert(cnt <= src_size);
			buckets[t][b] = cnt;
		}
	}
}

void
validate_ordered_src(u64* src, const u64 src_size)
{
	for (u64 i = 0; i != src_size - 1; ++i) {
		u64 hi = EXTRACT_HASH(src[i]);
		u64 hi1 = EXTRACT_HASH(src[i + 1]);
		r_assert(hi <= hi1)(i)(hi)(hi1);
	}
}

void
radix_sort(u64* src, const u64 src_size, const int num_threads)
{
	/// buckets
	u64** buckets = new u64*[num_threads];
	u64** tmp_buckets = new u64*[num_threads];
	for (int i = 0; i < num_threads; ++i) {
		buckets[i] = new u64[BS_BucketSize];
		fill(buckets[i], buckets[i] + BS_BucketSize, 0);
		tmp_buckets[i] = new u64[BS_BucketSize];
	}
	/// next buckets
	u64*** next_buckets = new u64**[num_threads];
	for (int t = 0; t < num_threads; ++t) {
		next_buckets[t] = new u64*[num_threads];
		for (int i = 0; i < num_threads; ++i) {
			next_buckets[t][i] = new u64[BS_BucketSize];
		}
	}
	
	BucketSortData* bsds = new BucketSortData[num_threads];
	int thread_id;
	pthread_mutex_t tid_lock;
	pthread_mutex_init(&tid_lock, NULL);
	pthread_t tids[num_threads];
	u64* trg = new u64[src_size];
	char process_msg[1024];
	
	for (u64 p = 0; p != BS_PASS; ++p) {		
		sprintf(process_msg, "sort pass %u", (unsigned)p);
		DynamicTimer dt(process_msg);
		if (p == 0) {
			init_buckets(buckets, num_threads, src, src_size);
			//init_buckets_mt(buckets, num_threads, src, src_size);
		} else {
			fill_bucket_from_next(buckets, next_buckets, num_threads);
		}
		
		update_bucket(buckets, tmp_buckets, src_size, num_threads);
		for (int t1 = 0; t1 < num_threads; ++t1)
			for (int t2 = 0; t2 < num_threads; ++t2)
				for (u64 b = 0; b != BS_BucketSize; ++b) {
					next_buckets[t1][t2][b] = 0;
				}
		
		thread_id = 0;
		for (int i = 0; i < num_threads; ++i) {
			bsds[i].src          = src;
			bsds[i].trg          = trg;
			bsds[i].src_size     = src_size;
			bsds[i].pass         = &p;
			bsds[i].gtid         = &thread_id;
			bsds[i].gtid_lock    = &tid_lock;
			bsds[i].num_threads  = num_threads;
			bsds[i].buckets      = buckets;
			bsds[i].next_buckets = next_buckets;
			
			pthread_create(tids + i, NULL, sort_thread, bsds + i);
		}
		for (int i = 0; i < num_threads; ++i) {
			pthread_join(tids[i], NULL);
		}
		
		u64* tmp = src;
		src = trg;
		trg = tmp;
	}
	
	for (int i = 0; i < num_threads; ++i) {
		delete[] buckets[i];
		delete[] tmp_buckets[i];
	}
	delete[] buckets;
	delete[] tmp_buckets;
	
	for (int i = 0; i < num_threads; ++i) {
		for (int j = 0; j < num_threads; ++j) {
			delete[] next_buckets[i][j];
		}
		delete[] next_buckets[i];
	}
	delete[] next_buckets;
	
	delete[] bsds;
	delete[] trg;
	validate_ordered_src(src, src_size);
}

void
build_kmer_starts(u64* offset_list, const u64 num_kmers, const int kmer_size, u64**& kmer_starts)
{
	DynamicTimer dt(__func__);
	
	u64 max_hash = 1; max_hash = max_hash << (kmer_size * 2);
	kmer_starts = new u64*[max_hash];
	fill(kmer_starts, kmer_starts + max_hash, (u64*)NULL);
	u64 i = 0, j;
	while (i < num_kmers) {
		const u64 hash = EXTRACT_HASH(offset_list[i]);
		j = i + 1;
		while (j < num_kmers) {
			const u64 h = EXTRACT_HASH(offset_list[j]);
			if (h != hash) break;
			++j;
		}
		kmer_starts[hash] = offset_list + i;
		i = j;
	}
}

void
clear_hash_in_offset_list(u64* offset_list, const u64 num_kmers, const u64 max_offset)
{
	DynamicTimer dt(__func__);
	
	for (u64 i = 0; i != num_kmers; ++i) {
		offset_list[i] = offset_list[i] & OFFSET_MASK;
		r_assert(offset_list[i] < max_offset);
	}
}

void 
build_lookup_table(const PackedDB* ref, 
				   const int kmer_size,
				   const int kmer_cnt_cutoff,
				   const int num_threads,
				   short*& kmer_cnts, 
				   u64**& kmer_starts, 
				   u64*& offset_list)
{
	u64 num_kmers;
	kmer_cnts = get_kmer_counts(ref, kmer_size, kmer_cnt_cutoff, num_kmers);
	offset_list = get_offset_list(ref, kmer_cnts, kmer_size, num_kmers);
	radix_sort(offset_list, num_kmers, num_threads);
	build_kmer_starts(offset_list, num_kmers, kmer_size, kmer_starts);
	clear_hash_in_offset_list(offset_list, num_kmers, ref->size());
}

} // ns_mc_lookup_table
