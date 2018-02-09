#include "consensus_one_partition.h"

#include "consensus_aux.h"
#include "consensus_one_read.h"
#include "../common/ontcns_aux.h"
#include "../partition_candidates/pcan_aux.h"

/*
static void
load_partition_candidates(const char* can_path, const int pid, int* min_read_id, int* max_read_id, vec_can* candidates)
{
	new_kstring(pname);
	make_partition_name(can_path, pid, &pname);
	char line[2048];
	kv_clear(*candidates);
	GappedCandidate can;
	int Lid = 1000000000, Rid = -1;
	DGZ_OPEN(in, kstr_str(pname), "r");
	while (gzgets(in, line, 2048)) {
		LOAD_GAPPED_CANDIDATE(sscanf, line, can);
		kv_push(GappedCandidate, *candidates, can);
		Lid = OC_MIN(Lid, can.sid);
		Rid = OC_MAX(Rid, can.sid);
	}
	GZ_CLOSE(in);
	free_kstring(pname);
	ks_introsort_GappedCandidate_SidLT(kv_size(*candidates), kv_data(*candidates));
	
	*min_read_id = Lid;
	*max_read_id = Rid + 1;
}
*/

static void
load_partition_candidates(const char* can_path, const int pid, int* min_read_id, int* max_read_id, vec_can* candidates)
{
	new_kstring(pname);
	make_partition_name(can_path, pid, &pname);
	size_t file_bytes = FILE_SIZE(kstr_data(pname));
	size_t num_can = file_bytes / sizeof(GappedCandidate);
	
	kv_clear(*candidates);
	kv_resize(GappedCandidate, *candidates, num_can);
	DFOPEN(in, kstr_data(pname), "rb");
	FREAD(kv_data(*candidates), sizeof(GappedCandidate), num_can, in);
	FCLOSE(in);
	free_kstring(pname);
	ks_introsort_GappedCandidate_SidLT(kv_size(*candidates), kv_data(*candidates));
	*min_read_id = kv_A(*candidates, 0).sid;
	*max_read_id = kv_A(*candidates, num_can - 1).sid + 1;
}

static void*
consensus_thread(void* arg)
{
	CnsData* cns_data = (CnsData*)(arg);
	size_t can_sid;
	size_t can_eid;
	while (extract_candidate_range(cns_data, &can_sid, &can_eid)) {
		consensus_one_read(cns_data, can_sid, can_eid);
		//if (can_sid > 100) break;
	}
	return NULL;
}

void
consensus_one_partition(PackedDB* reads,
						const char* can_path,
						CnsOptions* options,
						FILE* cns_out,
						FILE* raw_out,
						const int pid)
{
	char job_name[1024];
	sprintf(job_name, "consensus partition %d", pid);
	TIMING_START(job_name);
	
	size_t next_can_id = 0;
	OcMutex can_id_lock;
	pthread_mutex_init(&can_id_lock, NULL);
	OcMutex cns_out_lock;
	pthread_mutex_init(&cns_out_lock, NULL);
	OcMutex raw_out_lock;
	pthread_mutex_init(&raw_out_lock, NULL);
	OcMutex read_id_lock;
	pthread_mutex_init(&read_id_lock, NULL);
	ReadIdPool* corrected_read_ids = new_ReadIdPool(&read_id_lock);
	new_kvec(vec_can, candidates);
	int max_read_id, min_read_id;
	load_partition_candidates(can_path, pid, &min_read_id, &max_read_id, &candidates);
	
	CnsData** cns_data_pool = (CnsData**)malloc(options->num_threads * sizeof(CnsData*));
	for (int i = 0; i < options->num_threads; ++i) {
		cns_data_pool[i] = new_CnsData(reads,
									&candidates,
									&next_can_id,
									&can_id_lock,
									options,
									cns_out,
									&cns_out_lock,
									raw_out,
									&raw_out_lock,
									corrected_read_ids);
	}
	
	pthread_t* tids = (pthread_t*)malloc(sizeof(pthread_t) * options->num_threads);
	for (int i = 0; i < options->num_threads; ++i) {
		pthread_create(&tids[i], NULL, consensus_thread, cns_data_pool[i]);
	}
	for (int i = 0; i < options->num_threads; ++i) {
		pthread_join(tids[i], NULL);
	}
	free(tids);
	
	for (int i = 0; i < options->num_threads; ++i) {
		cns_data_pool[i] = free_CnsData(cns_data_pool[i]);
	}
	free(cns_data_pool);
	
	new_CnsSeq(cns_seq);
	for (int i = min_read_id; i < max_read_id; ++i) {
		if (!read_id_exists(corrected_read_ids, i)) {
			cns_seq.id = i;
			cns_seq.left = 0;
			cns_seq.right = PDB_SEQ_SIZE(reads, i);
			cns_seq.org_seq_size = cns_seq.right;
			pdb_extract_sequence(reads, i, FWD, &cns_seq.cns_seq);
			for (size_t k = 0; k != kstr_size(cns_seq.cns_seq); ++k) {
				kstr_A(cns_seq.cns_seq, k) = DecodeDNA(kstr_A(cns_seq.cns_seq, k));
			}
			kputc('\0', &cns_seq.cns_seq);
			DUMP_CNS_SEQ(fprintf, raw_out, cns_seq);
		}
	}
	free_CnsSeq(cns_seq);
	
	free_kvec(candidates);
	free_ReadIdPool(corrected_read_ids);
	
	TIMING_END(job_name);
}
