#include "read_id_pool.h"

#include "../klib/khash.h"

KHASH_SET_INIT_INT(32)
typedef khash_t(32) khash_int;

ReadIdPool*
new_ReadIdPool(OcMutex* id_lock)
{
	ReadIdPool* pool = (ReadIdPool*)malloc(sizeof(ReadIdPool));
	pool->read_id_pool = (void*)kh_init(32);
	pool->id_lock = id_lock;
	return pool;
}

ReadIdPool*
free_ReadIdPool(ReadIdPool* pool)
{
	kh_destroy(32, pool->read_id_pool);
	free(pool);
	return 0;
}

void
clear_ReadIdPool(ReadIdPool* pool)
{
	kh_clear_32(pool->read_id_pool);
}

BOOL
read_id_exists(ReadIdPool* pool, int id)
{
	khash_int* hid = (khash_int*)pool->read_id_pool;
	return kh_get(32, hid, id) != kh_end(hid);
}

void
add_read_id(ReadIdPool* pool, int id)
{
	khash_int* hid = (khash_int*)pool->read_id_pool;
	int ret;
	if (pool->id_lock) pthread_mutex_lock(pool->id_lock);
	kh_put(32, hid, id, &ret);
	if (pool->id_lock) pthread_mutex_unlock(pool->id_lock);
}
