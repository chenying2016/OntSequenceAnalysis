#ifndef OPTIONS_H
#define OPTIONS_H

#include "../word_finder/word_finder_aux.h"

struct MapOptions
{
    const char* wrk_dir;
    const char* output;
    int         num_threads;
    int         kmer_size;
	int			kmer_cnt_cutoff;
    int         scan_stride;
    int         block_size;
    int         num_candidates;
    int         block_score_cutoff;
	double      error_rate;
	
	InitWordParameters get_word_params() const {
		InitWordParameters p;
		p.kmer_size = kmer_size;
		p.scan_stride = scan_stride;
		p.num_candidates = num_candidates;
		p.block_score_cutoff = block_score_cutoff;
		p.block_size = block_size;
		p.error_rate = error_rate;
		p.ddfs_cutoff = error_rate;
		
		return p;
	}
};

int
parse_arguments(int argc, char* argv[], MapOptions& options);

void
print_usage();

void
print_options(const MapOptions& p);

#endif // OPTIONS_H
