#include "read_id_pool.h"

#include "cns_options.h"
#include "../common/makedb_aux.h"
#include "../partition_candidates/pcan_aux.h"
#include "consensus_one_partition.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	FILE* out = stderr;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s [options] wrk_dir candidates cns_out raw_out\n", prog);
	fprintf(out, "\n");
	fprintf(out, "OPTIONS AND DESCRIPTIONS:\n");
	describe_CnsOptions();
}

int main(int argc, char* argv[])
{
	CnsOptions options;
	if (argc < 5 || parse_CnsOptions(argc - 4, argv, &options) != ARG_PARSE_SUCCESS) {
		print_usage(argv[0]);
		return 1;
	}
	
	print_CnsOptions(&options);
	const char* wrk_dir = argv[argc - 4];
	const char* candidates_path = argv[argc - 3];
	const char* cns_out_path = argv[argc - 2];
	const char* raw_out_path = argv[argc - 1];
	PackedDB* reads = merge_volumes(wrk_dir);
	pdb_print_info(reads);
	int num_partitions = load_num_partitions(candidates_path);
	DFOPEN(cns_out, cns_out_path, "w");
	DFOPEN(raw_out, raw_out_path, "w");
	
	for (int i = 0; i < num_partitions; ++i) {
		consensus_one_partition(reads, candidates_path, &options, cns_out, raw_out, i);
	}
	
	FCLOSE(cns_out);
	FCLOSE(raw_out);
	free_PackedDB(reads);
}
