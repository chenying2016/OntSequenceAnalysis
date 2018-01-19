#include "results_writer.h"

ResultsWriter*
new_results_writer(FILE* _out, pthread_mutex_t* _out_lock)
{
	ResultsWriter* writer = (ResultsWriter*)malloc(sizeof(ResultsWriter));
	writer->buffer = new_kstring();
	writer->out = _out;
	writer->out_lock = _out_lock;
	
	return writer;
}

void
dump_results_writer(ResultsWriter* writer)
{
	size_t n = KSTR_SIZE(*writer->buffer);
	if (!n) return;
	if (writer->out_lock) {
		pthread_mutex_lock(writer->out_lock);
	}
	FWRITE(KSTR_DATA(*writer->buffer), 1, n, writer->out);
	if (writer->out_lock) {
		pthread_mutex_unlock(writer->out_lock);
	}
	KSTR_CLEAR(*writer->buffer);
}

void
check_results_writer(ResultsWriter* writer)
{
	const size_t kMaxBufferSize = 8 * (1 << 20);
	size_t n = KSTR_SIZE(*writer->buffer);
	if (n >= kMaxBufferSize) {
		dump_results_writer(writer);
	}
}

ResultsWriter* 
destroy_results_writer(ResultsWriter* writer)
{
	if (!writer) return 0;
	dump_results_writer(writer);
	free_kstring(writer->buffer);
	free(writer);
	return 0;
}
