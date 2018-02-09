#include "pcan.h"

#include <assert.h>

#include "pcan_aux.h"
#include "pcan_options.h"
#include "../common/makedb_aux.h"
#include "../common/ontcns_aux.h"
#include "../common/record_reader.h"

#define fix_candidate_offsets(c, nc, query_is_target) \
do { \
	nc.score = c.score; \
	\
	if (query_is_target) { \
		nc.qid = c.sid; \
		nc.qsize = c.ssize; \
		nc.sid = c.qid; \
		nc.ssize = c.qsize; \
		\
		if (c.qdir == REV) { \
			nc.qdir = ReverseStrand(c.sdir); \
			nc.qbeg = c.ssize - c.send; \
			nc.qend = c.ssize - c.sbeg; \
			nc.qoff = c.ssize - c.soff; \
			nc.sdir = FWD; \
			nc.sbeg = c.qsize - c.qend; \
			nc.send = c.qsize - c.qbeg; \
			nc.soff = c.qsize - c.qoff; \
		} else { \
			assert(c.qdir == FWD); \
			nc.qdir = c.sdir; \
			nc.qbeg = c.sbeg; \
			nc.qend = c.send; \
			nc.qoff = c.soff; \
			nc.sdir = c.qdir; \
			nc.sbeg = c.qbeg; \
			nc.send = c.qend; \
			nc.soff = c.qoff; \
		} \
	} else { \
		nc.qid = c.qid; \
		nc.qsize = c.qsize; \
		nc.sid = c.sid; \
		nc.ssize = c.ssize;	\
		\
		if (c.sdir == REV) {	\
			nc.qdir = ReverseStrand(c.qdir);	\
			nc.qbeg = c.qsize - c.qend; \
			nc.qend = c.qsize - c.qbeg; \
			nc.qoff = c.qsize - c.qoff; \
			nc.sdir = FWD; \
			nc.sbeg = c.ssize - c.send; \
			nc.send = c.ssize - c.sbeg; \
			nc.soff = c.ssize - c.soff; \
		} else { \
			assert(c.sdir == FWD); \
			nc = c; \
		} \
	} \
} while(0)

void
pcan_main(PcanOptions* options,
		  const char* wrk_dir,
		  const char* can_path)
{
	int num_reads = load_num_reads(wrk_dir);
	int num_batches = (num_reads + options->batch_size - 1) / options->batch_size;
	dump_num_partitions(can_path, num_batches);
	CandidatePartitionWriter* writer = new_CandidatePartitionWriter(options->num_output_files);
	char line[1024];
	GappedCandidate can, ncan;
	char job[1024];
	
	for (int fid = 0; fid < num_batches; fid += options->num_output_files) {
		int sfid = fid;
		int efid = OC_MIN(sfid + options->num_output_files, num_batches);
		int nfid = efid - sfid;
		int ssid = sfid * options->batch_size;
		int esid = efid * options->batch_size;
		sprintf(job, "dumping candidates for partitions [%d, %d)", sfid, efid);
		TIMING_START(job);
		open_CandidatePartitionWriter(nfid, can_path, sfid, writer);
		BinRecordReader* m4in = new_BinRecordReader(can_path, sizeof(GappedCandidate));
		while (get_BinRecordReader(m4in, &can)) {
			if (can.qid >= ssid && can.qid < esid) {
				fix_candidate_offsets(can, ncan, TRUE);
				int i = (can.qid - ssid) / options->batch_size;
				dump_GappedCandidate(writer, i, &ncan);
			}
			if (can.sid >= ssid && can.sid < esid) {
				fix_candidate_offsets(can, ncan, FALSE);
				int i = (can.sid - ssid) / options->batch_size;
				dump_GappedCandidate(writer, i, &ncan);
			}
		}
		close_CandidatePartitionWriter(writer);
		free_BinRecordReader(m4in);
		TIMING_END(job);
	}
	
	free_CandidatePartitionWriter(writer);
}
