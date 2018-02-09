#ifndef PCAN_AUX_H
#define PCAN_AUX_H

#include "../klib/kstring.h"
#include "../common/gapped_candidate.h"
#include "../common/record_writer.h"

void
make_partition_name(const char* prefix, const int pid, kstring_t* name);

void
make_partition_index_name(const char* prefix, kstring_t* name);

int 
load_num_partitions(const char* candidates_path);

void
dump_num_partitions(const char* candidates_path, const int np);

typedef struct {
	int max_output_files;
	int num_output_files;
	ResultsWriter** rlist;
} CandidatePartitionWriter;

CandidatePartitionWriter*
new_CandidatePartitionWriter(const int max_output_files);

CandidatePartitionWriter*
free_CandidatePartitionWriter(CandidatePartitionWriter* writer);

void
close_CandidatePartitionWriter(CandidatePartitionWriter* writer);

void
open_CandidatePartitionWriter(const int num_output_files, const char* prefix, const int sfid, CandidatePartitionWriter* writer);

void
dump_GappedCandidate(CandidatePartitionWriter* writer, const int fid, GappedCandidate* can);

#endif // PCAN_AUX_H
