#include "mc_consensus.h"

#include "consensus_aux.h"
#include "overlaps_partition.h"
#include "../common/cns_seq.h"
#include "../common/pdb_aux.h"
#include "../common/record_reader.h"
#include "../common/timer.h"
#include "cns_stage1.h"
#include "cns_stage2.h"

#include <algorithm>
#include <fstream>
#include <vector>

using namespace std;

void*
consensus_func(void* arg)
{
	ConsensusData& cns_data = *(ConsensusData*)(arg);
	size_t cani_sid, cani_eid;
	while (cns_data.next_cani_range(cani_sid, cani_eid)) {
		cns_stage2(&cns_data, cani_sid, cani_eid);
	}
	
	return NULL;
}

static void
load_candidates(const char* path, vector<GappedCandidate>& canv)
{
	canv.clear();
	RecordReader can_reader(path);
	GappedCandidate can;
	while (can_reader.read_one_record(can)) {
		canv.push_back(can);
	}
	sort(canv.begin(), canv.end(), GappedCandidateSid_LT());
}

void
consensus_one_partition(const char* input,
						const int min_tid,
						const int max_tid,
						CnsOptions& options,
						PackedDB& reads,
						std::ostream& out,
					    std::ostream& raw_out)
{
	vector<GappedCandidate> canv;
	load_candidates(input, canv);
	const int num_threads = options.num_threads;
	pthread_mutex_t out_lock;
	pthread_mutex_init(&out_lock, NULL);
	pthread_mutex_t raw_out_lock;
	pthread_mutex_init(&raw_out_lock, NULL);
	pthread_mutex_t can_id_lock;
	pthread_mutex_init(&can_id_lock, NULL);
	size_t next_can_idx = 0;
	pthread_t thread_ids[num_threads];
	ConsensusData*  pcds[num_threads];

	build_cns_thrd_data(&reads,
						&canv,
						&next_can_idx,
						&can_id_lock,
						&options,
						&out,
						&out_lock,
						&raw_out,
						&raw_out_lock,
						pcds);
						
	for (int i = 0; i < num_threads; ++i) {
		pthread_create(&thread_ids[i], NULL, consensus_func, (void*)(pcds[i]));
	}
	for (int i = 0; i < num_threads; ++i) {
		pthread_join(thread_ids[i], NULL);
	}
	
	///*
	set<int> corrected_read_ids;
	for (int i = 0; i < num_threads; ++i) {
		corrected_read_ids.insert(pcds[i]->corrected_read_ids.begin(), pcds[i]->corrected_read_ids.end());
	}
	vector<char> read;
	int from_id = min_tid;
	int to_id = min<idx>(max_tid, reads.num_seqs());
	CnsSeq cns_seq;
	for (int i = from_id; i < to_id; ++i) {
		if (corrected_read_ids.find(i) == corrected_read_ids.end()) {
			int size = reads.seq_size(i);
			read.resize(size + 1);
			reads.get_raw_sequence(i, 0, size, 1, read.data());
			
			cns_seq.id = extract_read_id_from_ontcns_hdr(reads.seq_header(i));
			cns_seq.left = 0;
			cns_seq.right = size;
			cns_seq.org_seq_size = size;
			cns_seq.seq.assign(read.begin(), read.begin() + size);
			raw_out << cns_seq;
		}
	}
	//*/
	
	destroy_cns_thrd_data(pcds, num_threads);
}

void
cns_main(CnsOptions& options)
{
	const int batch_size = options.batch_size;
	const int num_batches = partition_candidates(options); 
	PackedDB reads;
	load_packed_reads(options.wrk_dir, reads);
	string pname;
	dsopen(ofstream, out, options.corrected_reads, ios::out);
	string raw_out_name = options.uncorrected_reads;
	dsopen(ofstream, raw_out, raw_out_name.c_str(), ios::out);
	char process_msg[256];
	for (int i = 0; i < num_batches; ++i) {
		int min_tid = i * batch_size;
		int max_tid = min_tid + batch_size;
		make_partition_file_name(options.candidates, i, pname);
		sprintf(process_msg, "consensus %s", pname.c_str());
		DynamicTimer dt(process_msg);
		consensus_one_partition(pname.c_str(), min_tid, max_tid, options, reads, out, raw_out);
	}
	sclose(out);
	sclose(raw_out);
}
