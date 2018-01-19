#include "map_aux.h"

MappingThreadData*
new_mapping_thread_data(int thread_id,
						MapOptions* options,
						PackedDB* reads,
						int read_start_id,
						PackedDB* reference,
						int reference_start_id,
						LookupTable* lktbl,
						FILE* out_file,
						pthread_mutex_t* out_lock,
						const int chunk_size,
						int* chunk_id,
						pthread_mutex_t* chunk_lock)
{
	MappingThreadData* mtd = (MappingThreadData*)malloc(sizeof(MappingThreadData));
	mtd->thread_id			= thread_id;
	mtd->options 			= options;
	mtd->reads 				= reads;
	mtd->read_start_id 		= read_start_id;
	mtd->reference 			= reference;
	mtd->reference_start_id = reference_start_id;
	mtd->lktbl 				= lktbl;
	mtd->output     		= new_results_writer(out_file, out_lock);
	mtd->chunk_id 			= chunk_id;
	mtd->chunk_size 		= chunk_size;
	mtd->chunk_lock 		= chunk_lock;
	return mtd;
}

MappingThreadData*
destroy_mapping_thread_data(MappingThreadData* mtd)
{
	mtd->output = destroy_results_writer(mtd->output);
	free(mtd);
	return 0;
}

BOOL
get_next_read_chunk(const int chunk_size,
					int* chunk_id,
					pthread_mutex_t* chunk_lock,
					const int num_reads,
					int* sid,
					int* eid)
{
	int next_chunk_id;
	pthread_mutex_lock(chunk_lock);
	next_chunk_id = *chunk_id;
	++(*chunk_id);
	pthread_mutex_unlock(chunk_lock);
	int L = next_chunk_id * chunk_size;
	int R = L + chunk_size;
	if (L >= num_reads) return FALSE;
	*sid = L;
	*eid = (R > num_reads) ? num_reads : R;
	return TRUE;
}
