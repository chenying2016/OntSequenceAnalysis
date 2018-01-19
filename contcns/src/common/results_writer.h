#ifndef RESULTS_WRITER_H
#define RESULTS_WRITER_H

#include <pthread.h>

#include "gapped_candidate.h"
#include "kstring.h"
#include "ontcns_aux.h"

typedef struct
{
	kstring_t* 			buffer;
	FILE*				out;
	pthread_mutex_t*	out_lock;
} ResultsWriter;

ResultsWriter*
new_results_writer(FILE* _out, pthread_mutex_t* _out_lock);

void
dump_results_writer(ResultsWriter* writer);

void
check_results_writer(ResultsWriter* writer);

ResultsWriter* 
destroy_results_writer(ResultsWriter* writer);

#define DUMP_RESULT(writer, fmt, args...) \
	do { \
		ksprintf(writer->buffer, fmt, ##args); \
		check_results_writer(writer); \
	} while(0)
   
#define DUMP_GAPPED_CANDIDATE(writer, can) \
	do { \
		kstring_output_GappedCandidate(ksprintf, writer->buffer, can); \
		check_results_writer(writer); \
	} while(0)

#endif // RESULTS_WRITER_H
