#include "pm_worker.h"

#include "../common/map_aux.h"
#include "../common/makedb_aux.h"
#include "../lookup_table/lookup_table.h"
#include "../word_finder/word_finder.h"

void*
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
	
	kstring_t* read = new_kstring();
	kstring_t* subject = new_kstring();
	WordFindData* wfdata = new_WordFindData(options->kmer_size, options->block_score_cutoff);
	vec_can candidates; kv_init(candidates);
	
	int sid, eid;
	while (get_next_read_chunk(mtd->chunk_size, mtd->chunk_id, mtd->chunk_lock, num_reads, &sid, &eid)) {
		//OC_LOG("mapping read %d --- %d\n", sid, eid);
		for (int i = sid; i < eid; ++i) {
			//OC_LOG("mapping read %d", i);
			kv_clear(candidates);
			pdb_extract_sequence(reads, i, FWD, read);
			find_candidates(KSTR_DATA(*read),
							KSTR_SIZE(*read),
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
			pdb_extract_sequence(reads, i, REV, read);
			find_candidates(KSTR_DATA(*read),
							KSTR_SIZE(*read),
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
			size_t n_can = kv_size(candidates);
			for (size_t i = 0; i != n_can; ++i) {
				kv_A(candidates, i).qid += read_start_id;
				kv_A(candidates, i).sid += reference_start_id;
			}
			
			if (n_can > options->num_output) {
				ks_introsort(GappedCandidate, n_can, kv_data(candidates));
			} 
			for (int k = 0; k < n_can && k < options->num_output; ++k) {
				GappedCandidate* can = &kv_A(candidates, k);
				DUMP_GAPPED_CANDIDATE(out, can);
			}
		}
	}
	
	free_kstring(read);
	free_kstring(subject);
	free_WordFindData(wfdata);
	kv_destroy(candidates);
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
		mdata[i] = new_mapping_thread_data(i,
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
		mdata[i] = destroy_mapping_thread_data(mdata[i]);
	}
	
	lktbl = destroy_lookup_table(lktbl);
	reference = free_PackedDB(reference);
	volumes = destroy_volumes_info(volumes);
	FCLOSE(out);
}
