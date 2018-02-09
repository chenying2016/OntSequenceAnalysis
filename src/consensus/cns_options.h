#ifndef CNS_OPTONS_H
#define CNS_OPTONS_H

#include "../common/ontcns_defs.h"

typedef struct {
	int align_size_cutoff;
	int	high_accuracy_coverage_cutoff;
	int low_accuracy_coverage_cutoff;
	int	min_size;
	int full_consensus;
	double error;
	int num_threads;
} CnsOptions;

void
describe_CnsOptions();

BOOL
parse_CnsOptions(int argc, char* argv[], CnsOptions* options);

void
print_CnsOptions(const CnsOptions* options);

#endif // CNS_OPTONS_H
