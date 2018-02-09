#include "pm_worker.h"

#include "../common/m4_record.h"
#include "../common/map_aux.h"
#include "../common/makedb_aux.h"
#include "../common/oc_assert.h"
#include "../lookup_table/lookup_table.h"
#include "../word_finder/word_finder.h"
#include "../gapped_align/oc_aligner.h"

RESULT_WRITER_IMPL(GappedCandidate, GAPPED_CANDIDATE, static) 
RESULT_WRITER_IMPL(M4Record, M4_RECORD, static) 

static void
extend_candidates(GappedCandidate* cans,
				  const int ncan,
				  OcAlignData* align_data,
				  const char* fwd_read,
				  const char* rev_read,
				  kstring_t* subject,
				  PackedDB* reads,
				  const int min_align_size,
				  vec_m4* m4list)
{
	ks_introsort_GappedCandidate_ScoreGT(ncan, cans);
	kv_clear(*m4list);
	M4Record m4;
	for (int i = 0; i < ncan; ++i) {
		if (check_candidate_contain(m4list, cans + i)) continue;
		//stream_output_GappedCandidate(fprintf, stdout, cans + i);
		const char* read = (cans[i].qdir == FWD) ? fwd_read : rev_read;
		pdb_extract_sequence(reads, cans[i].sid, FWD, subject);
		align_data->qid = cans[i].qid;
		align_data->tid = cans[i].sid;
		BOOL r = onc_align(read, cans[i].qoff, cans[i].qsize, kstr_str(*subject), cans[i].soff, cans[i].ssize, align_data, min_align_size);
		if (r) {
			m4.qid = cans[i].qid;
			m4.sid = cans[i].sid;
			m4.ident_perc = align_data->ident_perc;
			m4.vscore = cans[i].score;
			m4.qdir = cans[i].qdir;
			m4.qoff = align_data->qoff;
			m4.qend = align_data->qend;
			m4.qext = cans[i].qoff;
			m4.qsize = cans[i].qsize;
			m4.sdir = FWD;
			m4.soff = align_data->toff;
			m4.send = align_data->tend;
			m4.sext = cans[i].soff;
			m4.ssize = cans[i].ssize;
			//DUMP_M4_RECORD(fprintf, stdout, m4);
			kv_push(M4Record, *m4list, m4);
		}
	}
}

static void*
pm_search_one_volume(void* arg)
{
	MappingThreadData* mtd 	= (MappingThreadData*)(arg);
	ResultsWriter* out 		= mtd->output;
	PackedDB* reads 		= mtd->reads;
	const int num_reads		= PDB_NUM_SEQS(reads);
	PackedDB* reference 	= mtd->reference;
	int read_start_id		= mtd->read_start_id;
	int reference_start_id	= mtd->reference_start_id;
	LookupTable* lktbl		= mtd->lktbl;
	MapOptions* options		= mtd->options;
	
	new_kstring(fwd_read);
	new_kstring(rev_read);
	new_kstring(subject);
	WordFindData* wfdata = new_WordFindData(options->kmer_size, options->block_score_cutoff);
	new_kvec(vec_can, candidates);
	new_kvec(vec_m4, m4list);
	OcAlignData* align_data = new_OcAlignData(options->error);
	
	int sid, eid;
	while (get_next_read_chunk(mtd->chunk_size, mtd->chunk_id, mtd->chunk_lock, num_reads, &sid, &eid)) {
		//OC_LOG("mapping read %d --- %d", sid, eid);
		for (int i = sid; i < eid; ++i) {
			//OC_LOG("mapping read %d", i);
			kv_clear(candidates);
			pdb_extract_sequence(reads, i, FWD, &fwd_read);
			find_candidates(kstr_str(fwd_read),
							kstr_size(fwd_read),
							i,
							FWD,
							read_start_id,
							reference_start_id,
							TRUE,
							reference,
							lktbl,
							options,
							wfdata,
							&candidates);
			pdb_extract_sequence(reads, i, REV, &rev_read);
			find_candidates(kstr_str(rev_read),
							kstr_size(rev_read),
							i,
							REV,
							read_start_id,
							reference_start_id,
							TRUE,
							reference,
							lktbl,
							options,
							wfdata,
							&candidates);
			
			if (options->job == MAP_JOB_ALN) {
				ks_introsort(GappedCandidate_ScoreGT, kv_size(candidates), kv_data(candidates));
				if (kv_size(candidates) > options->num_candidates) kv_resize(GappedCandidate, candidates, options->num_candidates);
				//OC_LOG("number of candidates: %lu", kv_size(candidates));
				extend_candidates(kv_data(candidates),
								  kv_size(candidates),
								  align_data,
								  kstr_str(fwd_read),
								  kstr_str(rev_read),
								  &subject,
								  reference,
								  options->align_size_cutoff,
								  &m4list);
				for (size_t k = 0; k != kv_size(m4list); ++k) {
					kv_A(m4list, k).qid += read_start_id;
					kv_A(m4list, k).sid += reference_start_id;
					M4Record* m = &kv_A(m4list, k);
					RW_DUMP_ONE_DATA(M4_RECORD, out, m);
				}
			} else {
				oc_assert(options->job == MAP_JOB_CAN);
				for (size_t k = 0; k != kv_size(candidates); ++k) {
					kv_A(candidates, k).qid += read_start_id;
					kv_A(candidates, k).sid += reference_start_id;
				}
				if (kv_size(candidates) > options->num_output) {
					ks_introsort(GappedCandidate_ScoreGT, kv_size(candidates), kv_data(candidates));
					kv_resize(GappedCandidate, candidates, options->num_output);
				} 
				for (size_t k = 0; k < kv_size(candidates); ++k) {
					GappedCandidate* can = &kv_A(candidates, k);
					RW_DUMP_ONE_DATA(GAPPED_CANDIDATE, out, can) ;
				}
			}
		}
	}
	
	free_kstring(fwd_read);
	free_kstring(rev_read);
	free_kstring(subject);
	free_WordFindData(wfdata);
	kv_destroy(candidates);
	kv_destroy(m4list);
	free_OcAlignData(align_data);
	return NULL;
}

void
pm_main(MapOptions* options, const int vid, const char* wrk_dir, const char* output)
{
	print_MapOptions(options);
	VolumesInfo* volumes = load_volumes_info(wrk_dir);
	const char* reference_path = vi_volume_name(volumes, vid);
	PackedDB* reference = new_PackedDB();
	pdb_load(reference, reference_path, TECH_PACBIO);
	LookupTable* lktbl = build_lookup_table(reference, options->kmer_size, options->kmer_cnt_cutoff, options->num_threads);
	const int reference_start_id = kv_A(volumes->read_start_id, vid);
	DFOPEN(out, output, "w");
	pthread_mutex_t out_lock;
	pthread_mutex_init(&out_lock, NULL);
	pthread_mutex_t chunk_lock;
	pthread_mutex_init(&chunk_lock, NULL);
	int chunk_id;
	const int chunk_size = 500;
	pthread_t jobs[options->num_threads];
	MappingThreadData* mdata[options->num_threads];
	char job[1024];
	
	for (int i = 0; i < options->num_threads; ++i) {
		mdata[i] = new_MappingThreadData(i,
										   options, 
										   NULL, 
										   0,
										   reference,
										   reference_start_id,
										   lktbl,
										   out,
										   &out_lock,
										   chunk_size,
										   &chunk_id,
										   &chunk_lock);
	}
	
	for (int i = vid; i < volumes->num_volumes; ++i) {
		sprintf(job, "pairwise mapping v%d vs v%d", i, vid);
		TIMING_START(job);
		const char* reads_path = vi_volume_name(volumes, i);
		PackedDB* reads = new_PackedDB();
		pdb_load(reads, reads_path, TECH_PACBIO);
		int read_start_id = kv_A(volumes->read_start_id, i);
		chunk_id = 0;
		for (int k = 0; k < options->num_threads; ++k) {
			mdata[k]->reads = reads;
			mdata[k]->read_start_id = read_start_id;
			pthread_create(jobs + k, NULL, pm_search_one_volume, mdata[k]);
		}
		for (int k = 0; k < options->num_threads; ++k) {
			pthread_join(jobs[k], NULL);
		}
		reads = free_PackedDB(reads);
		TIMING_END(job);
	}
	
	for (int i = 0; i < options->num_threads; ++i) {
		mdata[i] = free_MappingThreadData(mdata[i]);
	}
	
	lktbl = destroy_lookup_table(lktbl);
	reference = free_PackedDB(reference);
	volumes = destroy_volumes_info(volumes);
	FCLOSE(out);
}
