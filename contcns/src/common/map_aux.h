#ifndef MAP_AUX_H
#define MAP_AUX_H

#include "map_options.h"
#include "packed_db.h"
#include "results_writer.h"
#include "../lookup_table/lookup_table.h"

typedef struct
{
	int					thread_id;
	MapOptions*			options;
	PackedDB*			reads;
	int					read_start_id;
	PackedDB*			reference;
	int					reference_start_id;
	LookupTable*		lktbl;
	ResultsWriter*		output;
	int					chunk_size;
	int*				chunk_id;
	pthread_mutex_t*	chunk_lock;
} MappingThreadData;

BOOL
get_next_read_chunk(const int chunk_size,
					int* chunk_id,
					pthread_mutex_t* chunk_lock,
					const int num_reads,
					int* sid,
					int* eid);

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
						pthread_mutex_t* chunk_lock);

MappingThreadData*
destroy_mapping_thread_data(MappingThreadData* mtd);

#endif // MAP_AUX_H
