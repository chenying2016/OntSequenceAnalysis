#ifndef RECORD_WRITER_H
#define RECORD_WRITER_H

#include <pthread.h>

#include "ontcns_aux.h"
#include "../klib/kstring.h"

#define OcWriterBufferDefaultSize (8*(1<<20))

typedef struct {
    FILE*       out;
    OcMutex*    out_lock;
    kstring_t   buffer;
    size_t      buffer_size;
} AscRecordWriter;

AscRecordWriter*
new_AscRecordWriter(FILE* out, OcMutex* out_lock, size_t buffer_size);

AscRecordWriter*
free_AscRecordWriter(AscRecordWriter* writer);

void
dump_AscRecordWriter(AscRecordWriter* writer);

void
check_AscRecordWriter(AscRecordWriter* writer);

#define ASC_RW_IMPL(structName, funcName, SCOPE) \
    SCOPE void asc_write_##funcName(AscRecordWriter* wptr, structName * rptr) { \
        DUMP_##funcName(ksprintf, &wptr->buffer, *rptr); \
        check_AscRecordWriter(wptr); \
    }

typedef struct {
    FILE*       out;
    OcMutex*    out_lock;
    char*       buffer;
	char*		curr;
    size_t      record_size;
    size_t      max_num_records;
    size_t      i;
} BinRecordWriter;

BinRecordWriter*
new_BinRecordWriter(FILE* out, OcMutex* out_lock, size_t record_size, size_t buffer_size);

BinRecordWriter*
free_BinRecordWriter(BinRecordWriter* writer);

void
dump_BinRecordWriter(BinRecordWriter* writer);

void
check_BinRecordWriter(BinRecordWriter* writer);

void
write_one_bin_result(BinRecordWriter* writer, void* ptr);

typedef struct {
    AscRecordWriter* asc_writer;
    BinRecordWriter* bin_writer;
    BOOL bin_result;
} ResultsWriter;

ResultsWriter*
new_ResultsWriter(FILE* out, OcMutex* out_lock, BOOL bin_result, size_t record_size, size_t buffer_size);

ResultsWriter*
free_ResultsWriter(ResultsWriter* writer);

void
dump_ResultsWriter(ResultsWriter* writer);

#define RESULT_WRITER_IMPL(structName, funcName, SCOPE) \
	ASC_RW_IMPL(structName, funcName, SCOPE) \
	SCOPE void rw_dump_##funcName(ResultsWriter* wptr, structName * rptr) { \
		if (wptr->bin_result) { \
			write_one_bin_result(wptr->bin_writer, (void*)rptr); \
		} else { \
			asc_write_##funcName(wptr->asc_writer, rptr); \
		} \
	}
   
#define RW_DUMP_ONE_DATA(funcName, wptr, rptr)  rw_dump_##funcName(wptr, rptr)

#endif // RECORD_WRITER_H
