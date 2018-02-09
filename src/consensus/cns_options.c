#include "cns_options.h"

#include <stdio.h>
#include <getopt.h>

#include "../common/ontcns_aux.h"

static const char* argn_list = "a:x:y:l:f:e:t:";

static const CnsOptions sDefaultCnsOptions = {
	500,  	// minimal align size
	10, 	// minimal coverage for high accuracy
	20, 	// minimal coverage for low accuracy
	500, 	// minimal corrected read size
	0,		// full consensus
	0.5, 	// sequencing error
	1 		// number of threads
};

void
print_CnsOptions(const CnsOptions* options)
{
	FILE* out = stderr;
	fprintf(out, "-a %d ", options->align_size_cutoff);
	fprintf(out, "-x %d ", options->high_accuracy_coverage_cutoff);
	fprintf(out, "-y %d ", options->low_accuracy_coverage_cutoff);
	fprintf(out, "-l %d ", options->min_size);
	fprintf(out, "-f %d ", options->full_consensus);
	fprintf(out, "-e %f ", options->error);
	fprintf(out, "-t %d ", options->num_threads);
	fprintf(out, "\n");
}

BOOL
parse_CnsOptions(int argc, char* argv[], CnsOptions* options)
{
	*options = sDefaultCnsOptions;
	int c;
	while ((c = getopt(argc, argv, argn_list)) != -1) {
		switch (c) {
			case 'a':
				options->align_size_cutoff = atoi(optarg);
				break;
			case 'x':
				options->high_accuracy_coverage_cutoff = atoi(optarg);
				break;
			case 'y':
				options->low_accuracy_coverage_cutoff = atoi(optarg);
				break;
			case 'l':
				options->min_size = atoi(optarg);
				break;
			case 'f':
				options->full_consensus = atoi(optarg);
				break;
			case 'e':
				options->error = atof(optarg);
				break;
			case 't':
				options->num_threads = atoi(optarg);
				break;
			case '?':
				fprintf(stderr, "invalid option '%c'\n", (char)c);
				return ARG_PARSE_FAIL;
				break;
			case ':':
				fprintf(stderr, "argument to option '%c' is not provided\n", (char)c);
				return ARG_PARSE_FAIL;
				break;
			default:
				break;
		}
	}
	return ARG_PARSE_SUCCESS;
}

void
describe_CnsOptions()
{
	FILE* out = stderr;
	fprintf(out, "-a <Integer>\talign length cutoff\n");
	fprintf(out, "-x <Integer>\tminimal coverage for sequencing error <= 15%%\n");
	fprintf(out, "-y <Integer>\tminimal coverage for sequencing error > 15%%\n");
	fprintf(out, "-l <Integer>\tminimal length of corrected reads.\n");
	fprintf(out, "-f <0 or 1>\tfull consensus or not: 1 = yes, 0 = no\n");
	fprintf(out, "-e <Real>\tsequencing error\n");
	fprintf(out, "-t <Integer>\tnumber of cpu threads\n");
	
	fprintf(out, "\n");
	fprintf(out, "DEFAULT OPTIONS:\n");
	print_CnsOptions(&sDefaultCnsOptions);
}
