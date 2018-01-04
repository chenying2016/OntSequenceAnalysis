#ifndef CNS_OPTIONS_H
#define CNS_OPTIONS_H

struct CnsOptions
{
	const char* wrk_dir;
    const char* candidates;
    const char* raw_reads;
    const char* corrected_reads;
	const char* uncorrected_reads;
	const char* reads_count_path;
    int         num_threads;
    int         batch_size;
    int         min_cov;
    int         min_size;
    int         min_align_size;
	double		min_mapping_ratio;
    double      error_rate;
    int         tech;
	int			full_cns;
	int 		align_method;
	int			normalise_gaps;
};

int 
parse_arguments(int argc, char* argv[], CnsOptions& options);

void
print_options(const CnsOptions& options);

void
print_usage();

#endif //  CNS_OPTIONS_H
