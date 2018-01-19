#ifndef MAP_OPTIONS_H
#define MAP_OPTIONS_H

#include "ontcns_defs.h"

typedef struct {
    int     kmer_size;
    int     scan_window;
    int     kmer_cnt_cutoff;
    int     block_size;
    int     block_score_cutoff;
    int     num_candidates;
    int     align_size_cutoff;
	double 	ddfs_cutoff;
    double  error;
    int     num_output;
    int     num_threads;
} MapOptions;

void
print_MapOptions(const MapOptions* p);

BOOL
parse_MapOptions(int argc, char* argv[], MapOptions* options);

void
describe_MapOptions();

#endif // MAP_OPTIONS_H
