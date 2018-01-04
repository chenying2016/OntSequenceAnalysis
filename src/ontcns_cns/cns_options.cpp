#include "cns_options.h"

#include <cstdlib>
#include <unistd.h>

#include <sstream>

#include "../common/defs.h"
#include "../common/mc_log.h"
#include "../gapped_align/ontcns_aligner.h"

using namespace std;

static const char* argn_list = "t:p:a:r:c:l:e:q:hf:s:m:n";

void
print_options(const CnsOptions& p)
{
    const char sep = ' ';

    cerr << "-t " << p.num_threads << sep
		 << "-p " << p.batch_size << sep
		 << "-a " << p.min_align_size << sep
		 << "-r " << p.min_mapping_ratio << sep
		 << "-c " << p.min_cov << sep
		 << "-l " << p.min_size << sep
		 << "-e " << p.error_rate << sep
		 << "-f " << p.full_cns << sep
		 << "-m " << p.align_method << sep
		 << "-n " << p.normalise_gaps
		 << endl;
	
	if (p.wrk_dir) cerr << "working directory: " << p.wrk_dir << endl;
    if (p.candidates) cerr << "candidates: " << p.candidates << endl;;
    if (p.raw_reads) cerr << "raw_reads: " << p.raw_reads << endl;
    if (p.corrected_reads)  cerr << "corrected_reads: " << p.corrected_reads << endl;
	if (p.uncorrected_reads) cerr << "uncorrected reads: " << p.uncorrected_reads << endl;
	if (p.reads_count_path) cerr << "reads_count_path: " << p.reads_count_path << endl;
}

static const CnsOptions
sNanoporeDefaultOptions = {
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	1,
	100000,
	6,
	500,
	400,
	0.0,
	0.2,
	TECH_NANOPORE,
	0,
	ALIGN_METHOD_EDLIB,
	0
};

template <typename T>
void parse_argument(const char opt_name, const char* arg, T& v)
{
    istringstream is(arg);
    if (!(is >> v)) mc_error << "failed to parse argument '" << arg << "' for option '" << opt_name << "'" << eolog;
}

/*
 * return value:
 * 0: parse successully
 * 1: only print usage
 * -1: parse error
 */

int
parse_arguments(int argc, char* argv[], CnsOptions& options)
{
    bool s_num_threads = 0;
    bool s_batch_size = 0;
    bool s_min_cov = 0;
    bool s_min_size = 0;
    bool s_min_align_size = 0;
	bool s_min_mapping_ratio = 0;
    bool s_error_rate = 0;
	bool s_min_ident = 0;
    bool s_tech = 0;
	bool s_full_cns = 0;
	bool s_stage = 0;
	bool print_usage = 0;
	bool s_align_method = 0;
	bool s_normalise_gaps = 0;
	int status = 0;
	
	CnsOptions tmp_opt;
	tmp_opt.tech = TECH_PACBIO;

    int c;
    while ((c = getopt(argc, argv, argn_list)) != -1) {
        switch (c) {
            case 't':
                parse_argument('t', optarg, tmp_opt.num_threads);
				if (tmp_opt.num_threads <= 0) {
					cerr << "number of threads (-t) must be > 0: " << tmp_opt.num_threads << endl;
					status = -1;
				}
				s_num_threads = 1;
                break;
            case 'p':
                parse_argument('k', optarg, tmp_opt.batch_size);
				if (tmp_opt.batch_size < 10000) {
					cerr << "batch size (-p) must be > 10000: " << tmp_opt.batch_size << endl;
					status = -1;
				}
				s_batch_size = 1;
                break;
            case 'a':
                parse_argument('a', optarg, tmp_opt.min_align_size);
				if (tmp_opt.min_align_size < 0) {
					cerr << "align size cutoff (-a) must be >= 0: " << tmp_opt.min_align_size << endl;
					status = -1;
				}
				s_min_align_size = 1;
                break;
			case 'r':
				parse_argument('r', optarg, tmp_opt.min_mapping_ratio);
				s_min_mapping_ratio = 1;
				break;
            case 'c':
                parse_argument('n', optarg, tmp_opt.min_cov);
				if (tmp_opt.min_cov < 1) {
					cerr << "coverage cutoff (-c) must be > 0: " << tmp_opt.min_cov << endl;
					status = -1;
				}
				s_min_cov = 1;
                break;
            case 'l':
                parse_argument('s', optarg, tmp_opt.min_size);
				if (tmp_opt.min_size < 1) {
					cerr << "corrected sequence length (-l) must be > 0: " << tmp_opt.min_size << endl;
					status = -1;
				}
				s_min_size = 1;
                break;
            case 'e':
                parse_argument('e', optarg, tmp_opt.error_rate);
				if (tmp_opt.error_rate < 0 || tmp_opt.error_rate >= 1.0) {
					cerr << "error rate (-e) must be in [0, 1.0]: " << tmp_opt.error_rate << endl;
					status = -1;
				}
				s_error_rate = 1;
                break;
			case 'f':
				parse_argument('f', optarg, tmp_opt.full_cns);
				if (tmp_opt.full_cns != 0 && tmp_opt.full_cns != 1) {
					cerr << "argument to option 'f' must be either '0' or '1': " << tmp_opt.full_cns << endl;
					status = -1;
				}
				s_full_cns = 1;
				break;
            case 'h':
				print_usage = 1;
                break;
			case 'm':
				parse_argument('m', optarg, tmp_opt.align_method);
				if (tmp_opt.align_method != ALIGN_METHOD_EDLIB
					&&
					tmp_opt.align_method != ALIGN_METHOD_DIFF) {
					cerr << "argument to option 'm' must be either '" << ALIGN_METHOD_EDLIB
						 << "' or '" << ALIGN_METHOD_DIFF << "'\n";
					status = -1;
				}
				s_align_method = 1;
				break;
			case 'n':
				tmp_opt.normalise_gaps = 1;
				s_normalise_gaps = 1;
				break;
			case '?':
				cerr << "invalid option '" << (char)c << "'" << endl;
				status = -1;
				break;
			case ':':
				cerr << "argument to option '" << (char)c << "' is not provided" << endl;
				status = -1;
            default:
                break;
        }
    }
	
	if (status == -1) return status;
	if (print_usage) return 1;
	if (argc - optind != 4) return -1;
	
	options = sNanoporeDefaultOptions;
	
	if (s_num_threads) options.num_threads = tmp_opt.num_threads;
    if (s_batch_size) options.batch_size = tmp_opt.batch_size;
    if (s_min_cov) options.min_cov = tmp_opt.min_cov;
    if (s_min_size) options.min_size = tmp_opt.min_size;
    if (s_min_align_size) options.min_align_size = tmp_opt.min_align_size;
	if (s_min_mapping_ratio) options.min_mapping_ratio = tmp_opt.min_mapping_ratio;
    if (s_error_rate) options.error_rate = tmp_opt.error_rate;
	if (s_full_cns) options.full_cns = tmp_opt.full_cns;
	if (s_align_method) options.align_method = tmp_opt.align_method;
	if (s_normalise_gaps) options.normalise_gaps = tmp_opt.normalise_gaps;
	
	options.wrk_dir = argv[optind];
	options.candidates = argv[optind + 1];
    //options.raw_reads = argv[optind + 1];
	options.corrected_reads = argv[optind + 2];
	options.uncorrected_reads = argv[optind + 3];
	
    return 0;
}

void
print_usage()
{
    cerr << "OntCns2Cns: a fast and accurate Nanopore sequence correction tool" << endl;
    cerr << endl << "USAGE: " << endl;
	cerr << "  " << "ontcns_cns" << " [options] " << "wrk-dir overlap-candidates corrected-reads uncorrected-reads" << endl;
	
	cerr << "\nOPTIONS:\n";
	cerr << "  -t <Integer>\t number of threads (cpus)" << endl;
	cerr << "  -b <Integer>\t batch size that reads will be partitioned" << endl;
	cerr << "  -a <Integer>\t align size cutoff" << endl;
	cerr << "  -c <Integer>\t coverage cutoff" << endl;
	cerr << "  -l <Integer>\t min lenght of corrected sequence" << endl;
	cerr << "  -e <Real>\t error rate of raw sequences" << endl;
	cerr << "  -f <0/1>\t full consensus (1) or not (0)" << endl;
	cerr << "  -s <1/2> consensus stage" << endl;
	cerr << "  -m <0/1> align method: 0 = edlib, 1 = diff" << endl;
	cerr << "  -n <0/1> normalise gaps: 0 = no, 1 = yes" << endl;
	cerr << "  -h \t\t print USAGE" << endl;
	
    cerr << "\nDefault options:" << endl << endl;
	print_options(sNanoporeDefaultOptions);
}
