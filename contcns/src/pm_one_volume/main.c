#include "pm_worker.h"
#include "../common/map_options.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	fprintf(stderr, "USAGE:\n");
	fprintf(stderr, "%s [options] wrk-dir volume-id output\n", prog);
	fprintf(stderr, "\n");
	fprintf(stderr, "OPTIONS AND DESCRIPTIONS:\n");
	describe_MapOptions();
}

int main(int argc, char* argv[])
{
	if (argc < 4) {
		print_usage(argv[0]);
		return 1;
	}
	MapOptions options;
	if (parse_MapOptions(argc - 3, argv, &options) != ARG_PARSE_SUCCESS) {
		print_usage(argv[0]);
		return 1;
	}
	const char* wrk_dir = argv[argc - 3];
	int vid = atoi(argv[argc - 2]);
	const char* output = argv[argc - 1];
	pm_main(&options, vid, wrk_dir, output);
	return 0;
}
