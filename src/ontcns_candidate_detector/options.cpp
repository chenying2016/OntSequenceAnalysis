#include "options.h"

#include <cstdlib>
#include <unistd.h>

#include <sstream>

#include "../common/defs.h"
#include "../common/mc_log.h"

using namespace std;

static const char* argn_list = "t:q:k:z:b:n:s:e:h";

void
print_options(const MapOptions& p)
{
    const char sep = ' ';

    cerr << "-t " << p.num_threads << sep
         << "-k " << p.kmer_size << sep
		 << "-q " << p.kmer_cnt_cutoff << sep
         << "-z " << p.scan_stride << sep
         << "-b " << p.block_size << sep
         << "-n " << p.num_candidates << sep
         << "-s " << p.block_score_cutoff << sep
         << "-e " << p.error_rate << sep
         << endl;
        
    if (p.wrk_dir) cerr << "wrk_dir: " << p.wrk_dir << endl;;
    if (p.output)  cerr << "output: " << p.output << endl;
}

static const MapOptions 
sNanoporeDefaultOptions_PM = {
    NULL, // wrk_dir
    NULL, // output
    1, // threads
    13, // kmer size
	1000, // kmer count cutoff
    10, // scan stride
    2000, // block size
    200, // number of candidates
    3, // block score cutoff
    0.5 // error 
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
parse_arguments(int argc, char* argv[], MapOptions& options)
{
    bool s_num_threads = 0;
    bool s_kmer_size = 0;
	bool s_kmer_cnt_cutoff = 0;
    bool s_scan_stride = 0;
    bool s_block_size = 0;
    bool s_num_candidates = 0;
    bool s_block_score_cutoff = 0;
    bool s_error_rate = 0;
	bool s_print_usage = 0;
	int status = 0;
	MapOptions tmp_opt;

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
            case 'k':
                parse_argument('k', optarg, tmp_opt.kmer_size);
				if (tmp_opt.kmer_size < 5) {
					cerr << "kmer size (-k) must be >= 5: " << tmp_opt.kmer_size << endl;
					status = -1;
				}
				s_kmer_size = 1;
                break;
			case 'q':
				parse_argument('q', optarg, tmp_opt.kmer_cnt_cutoff);
				s_kmer_cnt_cutoff = 1;
				break;
            case 'z':
                parse_argument('z', optarg, tmp_opt.scan_stride);
				if (tmp_opt.scan_stride < 0) {
					cerr << "scanning window (-z) must be >= 0: " << tmp_opt.scan_stride << endl;
					status = -1;
				}
				s_scan_stride = 1;
                break;
            case 'b':
                parse_argument('b', optarg, tmp_opt.block_size);
				if (tmp_opt.block_size < 500) {
					cerr << "block size (-b) must be >= 500: " << tmp_opt.block_size << endl;
					status = -1;
				}
				s_block_size = 1;
                break;
            case 'n':
                parse_argument('n', optarg, tmp_opt.num_candidates);
				if (tmp_opt.num_candidates <= 0) {
					cerr << "number of candidates (-n) must be > 0: " << tmp_opt.num_candidates << endl;
					status = -1;
				}
				s_num_candidates = 1;
                break;
            case 's':
                parse_argument('s', optarg, tmp_opt.block_score_cutoff);
				if (tmp_opt.block_score_cutoff <= 0) {
					cerr << "block score cutoff (-s) must be > 0: " << tmp_opt.block_score_cutoff << endl;
					status = -1;
				}
				s_block_score_cutoff = 1;
                break;
            case 'e':
                parse_argument('e', optarg, tmp_opt.error_rate);
				if (tmp_opt.error_rate < 0 || tmp_opt.error_rate >= 1.0) {
					cerr << "error rate (-e) must be in [0, 1.0]: " << tmp_opt.error_rate << endl;
					status = -1;
				}
				s_error_rate = 1;
                break;
            case 'h':
				s_print_usage = 1;
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
	
	if (status == -1) return -1;
	if (s_print_usage) return 1;
	if (argc - optind != 2) return -1;
	
	const char* wrk_dir = argv[optind];
	const char* output = argv[optind + 1];
	
	options = sNanoporeDefaultOptions_PM;
	
	if (s_num_threads) options.num_threads = tmp_opt.num_threads;
    if (s_kmer_size) options.kmer_size = tmp_opt.kmer_size;
	if (s_kmer_cnt_cutoff) options.kmer_cnt_cutoff = tmp_opt.kmer_cnt_cutoff;
    if (s_scan_stride) options.scan_stride = tmp_opt.scan_stride;
    if (s_block_size) options.block_size = tmp_opt.block_size;
    if (s_num_candidates) options.num_candidates = tmp_opt.num_candidates;
    if (s_block_score_cutoff) options.block_score_cutoff = tmp_opt.block_score_cutoff;
    if (s_error_rate) options.error_rate = tmp_opt.error_rate;
	
	options.wrk_dir = wrk_dir;
	options.output = output;

    return 0;
}

void
print_usage()
{
	cerr << "OntCns2CanFinder: a fast and accurate candidate detector for noisy long reads" << endl;
    cerr << endl << "USAGE: " << endl;
    cerr << "OntCns2CanFinder [options] wrk_dir output" << endl;

    cerr << "\nOPTIONS:\n";
    cerr << "  -t <Integer>\t number of threads (cpus)" << endl;
	cerr << "  -q <Integer>\t kmer occurs > q times will be ignored" << endl;
    cerr << "  -k <Integer>\t kmer size" << endl;
    cerr << "  -z <Integer>\t scanning window size" << endl;
    cerr << "  -b <Integer>\t block size" << endl;
    cerr << "  -n <Integer>\t number of candidate" << endl;
    cerr << "  -s <Integer>\t min number of kmer matches a matched block pair shares" << endl;
    cerr << "  -e <Real>\t sequencing error rate" << endl;
	cerr << "  -h\t\t print USAGE" << endl;

    cerr << "\nDefault options:" << endl << endl;
    print_options(sNanoporeDefaultOptions_PM);
}
