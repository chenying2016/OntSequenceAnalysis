#include "../common/map_options.h"
#include "../common/makedb_aux.h"
#include "../common/kstring.h"
#include "../common/ontcns_aux.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	fprintf(stderr, "USAGE:\n");
	fprintf(stderr, "%s [options] wrk-dir output\n", prog);
	fprintf(stderr, "\n");
	fprintf(stderr, "OPTIONS AND DESCRIPTIONS:\n");
	describe_MapOptions();
}

void
make_volume_result_name(const char* wrk_dir, const int vid, kstring_t* name)
{
	KSTR_CLEAR(*name);
	size_t n = strlen(wrk_dir);
	kputs(wrk_dir, name);
	if (wrk_dir[n - 1] != '/') kputc('/', name);
	ksprintf(name, "pm_result_%d", vid);
	kputc('\0', name);
}

void
make_pm_command(MapOptions* options, const char* wrk_dir, const int vid, const char* output, kstring_t* cmd)
{
	KSTR_CLEAR(*cmd);
	kputs("oc2pmov ", cmd);
	ksprintf(cmd, "-k %d -z %d -q %d -b %d -s %d -n %d -a %d -d %f -e %f -m %d -t %d ",
			 options->kmer_size,
			 options->scan_window,
			 options->kmer_cnt_cutoff,
			 options->block_size,
			 options->block_score_cutoff,
			 options->num_candidates,
			 options->align_size_cutoff,
			 options->ddfs_cutoff,
			 options->error,
			 options->num_output,
			 options->num_threads);
	ksprintf(cmd, "%s ", wrk_dir);
	ksprintf(cmd, "%d ", vid);
	ksprintf(cmd, "%s", output);
	kputc('\0', cmd);
}

void
make_merge_results_command(const int vid, const char* volume_results_name, const char* output, kstring_t* cmd)
{
	KSTR_CLEAR(*cmd);
	if (vid == 0) {
		ksprintf(cmd, "cat %s > %s", volume_results_name, output);
	} else {
		ksprintf(cmd, "cat %s >> %s", volume_results_name, output);
	}
	kputc('\0', cmd);
}

int main(int argc, char* argv[])
{
	TIMING_START(__func__);
	
	MapOptions options;
	if (argc < 3 || (parse_MapOptions(argc - 2, argv, &options) != ARG_PARSE_SUCCESS)) {
		print_usage(argv[0]);
		return 1;
	}
	
	const char* wrk_dir = argv[argc - 2];
	const char* output = argv[argc - 1];
	int num_volumes = load_num_volumes(wrk_dir);
	
	kstring_t* cmd = new_kstring();
	kstring_t* result_name = new_kstring();
	for (int i = 0; i < num_volumes; ++i) {
		make_volume_result_name(wrk_dir, i, result_name);
		make_pm_command(&options, wrk_dir, i, KSTR_DATA(*result_name), cmd);
		const char* cmd_str = KSTR_DATA(*cmd);
		fprintf(stdout, "Running command '%s'\n", cmd_str);
		SYSTEM(cmd_str);
	}
	
	for (int i = 0; i < num_volumes; ++i) {
		make_volume_result_name(wrk_dir, i, result_name);
		make_merge_results_command(i, KSTR_DATA(*result_name), output, cmd);
		const char* cmd_str = KSTR_DATA(*cmd);
		SYSTEM(cmd_str);
	}
	
	cmd = free_kstring(cmd);
	result_name = free_kstring(result_name);
	TIMING_END(__func__);
}
