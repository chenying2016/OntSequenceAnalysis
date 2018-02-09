#ifndef READ_ID_POOL_H
#define READ_ID_POOL_H

#include "../common/ontcns_defs.h"

typedef struct {
	void* read_id_pool;
	pthread_mutex_t* id_lock;
} ReadIdPool;

ReadIdPool*
new_ReadIdPool(OcMutex* id_lock);

ReadIdPool*
free_ReadIdPool(ReadIdPool* pool);

void
clear_ReadIdPool(ReadIdPool* pool);

BOOL
read_id_exists(ReadIdPool* pool, int id);

void
add_read_id(ReadIdPool* pool, int id);

#endif // READ_ID_POOL_H
