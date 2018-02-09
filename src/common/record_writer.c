#include "record_writer.h"

AscRecordWriter*
new_AscRecordWriter(FILE* out, OcMutex* out_lock, size_t buffer_size)
{
	AscRecordWriter* w = (AscRecordWriter*)malloc( sizeof(AscRecordWriter) );
	kstr_init(w->buffer);
	ks_reserve(&w->buffer, buffer_size * 1.5);
	w->out = out;
	w->out_lock = out_lock;
	return w;
}

AscRecordWriter*
free_AscRecordWriter(AscRecordWriter* writer)
{
	dump_AscRecordWriter(writer);
	free_kstring(writer->buffer);
	free(writer);
	return 0;
}

void
dump_AscRecordWriter(AscRecordWriter* writer)
{
	size_t n = kstr_size(writer->buffer);
	if (!n) return;
	if (writer->out_lock) {
		pthread_mutex_lock(writer->out_lock);
	}
	FWRITE(kstr_str(writer->buffer), 1, n, writer->out);
	if (writer->out_lock) {
		pthread_mutex_unlock(writer->out_lock);
	}
	kstr_clear(writer->buffer);
}

void
check_AscRecordWriter(AscRecordWriter* writer)
{
	size_t n = kstr_size(writer->buffer);
	if (n >= writer->buffer_size) {
		dump_AscRecordWriter(writer);
	}
}

BinRecordWriter*
new_BinRecordWriter(FILE* out, OcMutex* out_lock, size_t record_size, size_t buffer_size)
{
	BinRecordWriter* w = (BinRecordWriter*)malloc( sizeof(BinRecordWriter) );
	w->record_size = record_size;
	w->max_num_records = (buffer_size / record_size);
	buffer_size = w->max_num_records * record_size;
	w->buffer = (char*)malloc(buffer_size);
	
	w->out = out;
	w->out_lock = out_lock;
	w->curr = w->buffer;
	w->i = 0;
	
	return w;
}

BinRecordWriter*
free_BinRecordWriter(BinRecordWriter* writer)
{
	dump_BinRecordWriter(writer);
	free(writer->buffer);
	free(writer);
	return 0;
}

void
dump_BinRecordWriter(BinRecordWriter* writer)
{
	size_t n = writer->record_size * writer->i;
	if (!n) return;
	if (writer->out_lock) {
		pthread_mutex_lock(writer->out_lock);
	}
	FWRITE(writer->buffer, 1, n, writer->out);
	if (writer->out_lock) {
		pthread_mutex_unlock(writer->out_lock);
	}
	
	writer->curr = writer->buffer;
	writer->i = 0;
}

void
check_BinRecordWriter(BinRecordWriter* writer)
{
	if (writer->i >= writer->max_num_records) {
		dump_BinRecordWriter(writer);
	}
}

void
write_one_bin_result(BinRecordWriter* writer, void* ptr)
{
	check_BinRecordWriter(writer);
	memcpy(writer->curr, ptr, writer->record_size);
	++writer->i;
	writer->curr += writer->record_size;
}

ResultsWriter*
new_ResultsWriter(FILE* out, OcMutex* out_lock, BOOL bin_result, size_t record_size, size_t buffer_size)
{
	ResultsWriter* w = (ResultsWriter*)malloc( sizeof(ResultsWriter) );
	w->asc_writer = 0;
	w->bin_writer = 0;
	w->bin_result = bin_result;
	
	if (bin_result) {
		w->bin_writer = new_BinRecordWriter(out, out_lock, record_size, buffer_size);
	} else {
		w->asc_writer = new_AscRecordWriter(out, out_lock, buffer_size);
	}
	
	return w;
}

ResultsWriter*
free_ResultsWriter(ResultsWriter* writer)
{
	if (writer->bin_writer) {
		free_BinRecordWriter(writer->bin_writer);
	}
	
	if (writer->asc_writer) {
		free_AscRecordWriter(writer->asc_writer);
	}
	
	free(writer);
	return 0;
}

void
dump_ResultsWriter(ResultsWriter* writer)
{
	if (writer->asc_writer) {
		dump_AscRecordWriter(writer->asc_writer);
	}
	if (writer->bin_writer) {
		dump_BinRecordWriter(writer->bin_writer);
	}
}
