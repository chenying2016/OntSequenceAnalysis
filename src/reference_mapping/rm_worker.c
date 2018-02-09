#include "rm_worker.h"

#include "../common/m4_record.h"
#include "../common/map_aux.h"
#include "../common/makedb_aux.h"
#include "../common/oc_assert.h"
#include "../lookup_table/lookup_table.h"
#include "../word_finder/word_finder.h"
#include "../gapped_align/oc_aligner.h"

RESULT_WRITER_IMPL(M4Record, M4_RECORD, static)

PackedDB*
load_mapping_reads(kseq_t* read)
{
	if (kseq_read(read) < 0) return 0;
	const idx VolumeSize = 2000000000;
	
	PackedDB* reads = new_PackedDB();
	pdb_enlarge_size(reads, VolumeSize + 2000000);
	while (1) {
		pdb_add_one_seq(reads, read, TECH_NANOPORE);
		if (PDB_SIZE(reads) >= VolumeSize) break;
		if (kseq_read(read) < 0) break;
	}
	return reads;
}

static void
calc_reference_range(GappedCandidate* can, idx* from, idx* to, int* soff)
{
	idx n = (can->qoff < can->soff) ? (can->qoff * 1.3) : (can->soff);
	n = OC_MIN(n, can->soff);
	*from = can->soff - n;
	*soff = n;
	
	idx sr = can->ssize - can->soff;
	idx qr = can->qsize - can->qoff;
	n = (qr < sr) ? (qr * 1.3) : (sr);
	n = OC_MIN(n, sr);
	*to = can->soff + n;
}

static void
rm_extend_candidates(GappedCandidate* cans,
				  const int ncan,
				  OcAlignData* align_data,
				  const char* fwd_read,
				  const char* rev_read,
				  kstring_t* subject,
				  PackedDB* reads,
				  PackedDB* reference,
				  const int min_align_size,
				  vec_m4* m4list)
{
	ks_introsort_GappedCandidate_ScoreGT(ncan, cans);
	kv_clear(*m4list);
	M4Record m4;
	idx sfrom, sto;
	int soff, ssize;
	for (int i = 0; i < ncan; ++i) {
		if (check_candidate_contain(m4list, cans + i)) continue;
		const char* read = (cans[i].qdir == FWD) ? fwd_read : rev_read;
		calc_reference_range(cans + i, &sfrom, &sto, &soff);
		ssize = sto - sfrom;
		pdb_extract_subsequence(reference, cans[i].sid, sfrom, sto, FWD, subject);
		align_data->qid = cans[i].qid;
		align_data->tid = cans[i].sid;
		BOOL r = onc_align(read, cans[i].qoff, cans[i].qsize, kstr_str(*subject), soff, ssize, align_data, min_align_size);
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
			m4.soff = align_data->toff + sfrom;
			m4.send = align_data->tend + sfrom;
			m4.sext = cans[i].soff + sfrom;
			m4.ssize = cans[i].ssize;
			
			kv_push(M4Record, *m4list, m4);
		}
	}
}

static void*
rm_search_one_volume(void* arg)
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
							FALSE,
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
							FALSE,
							reference,
							lktbl,
							options,
							wfdata,
							&candidates);
			
			ks_introsort(GappedCandidate_ScoreGT, kv_size(candidates), kv_data(candidates));
			if (kv_size(candidates) > options->num_candidates) kv_resize(GappedCandidate, candidates, options->num_candidates);
				//OC_LOG("number of candidates: %lu", kv_size(candidates));
			rm_extend_candidates(kv_data(candidates),
								  kv_size(candidates),
								  align_data,
								  kstr_str(fwd_read),
								  kstr_str(rev_read),
								  &subject,
								  reads,
								  reference,
								  options->align_size_cutoff,
								  &m4list);
			for (size_t k = 0; k != kv_size(m4list); ++k) {
				kv_A(m4list, k).qid += read_start_id;
				kv_A(m4list, k).sid += reference_start_id;
				M4Record* m = &kv_A(m4list, k);
				RW_DUMP_ONE_DATA(M4_RECORD, out, m);
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
rm_main(MapOptions* options, const char* reads_path, const char* reference_path, const char* output)
{
	PackedDB* reference = new_PackedDB();
	pdb_load(reference, reference_path, TECH_NANOPORE);
	LookupTable* lktbl = build_lookup_table(reference, options->kmer_size, options->kmer_cnt_cutoff, options->num_threads);
	const int reference_start_id = 0;
	
	OcMutex out_lock;
	pthread_mutex_init(&out_lock, NULL);
	OcMutex chunk_lock;
	pthread_mutex_init(&chunk_lock, NULL);
	const int num_threads = options->num_threads;
	pthread_t tids[num_threads];
	int read_start_id = 0;
	int chunk_id;
	DFOPEN(out, output, "w");
	MappingThreadData* mdata[num_threads];
	for (int i = 0; i < num_threads; ++i) {
		mdata[i] = new_MappingThreadData(i,
										   options,
										   NULL,
										   read_start_id,
										   reference,
										   reference_start_id,
										   lktbl,
										   out,
										   &out_lock,
										   500,
										   &chunk_id,
										   &chunk_lock);
	}
	
	DGZ_OPEN(reads_in, reads_path, "r");
	kseq_t* read = kseq_init(reads_in);
	char job[1024];
	while (1) {
		PackedDB* reads = load_mapping_reads(read);
		if (reads == NULL) break;
		sprintf(job, "mapping %lu reads", PDB_NUM_SEQS(reads));
		TIMING_START(job);
		chunk_id = 0;
		for (int i = 0; i < num_threads; ++i) {
			mdata[i]->reads = reads;
			mdata[i]->read_start_id = read_start_id;
			pthread_create(tids + i, NULL, rm_search_one_volume, mdata[i]);
		}
		for (int i = 0; i < num_threads; ++i) {
			pthread_join(tids[i], NULL);
		}
		read_start_id += PDB_NUM_SEQS(reads);
		free_PackedDB(reads);
		TIMING_END(job);
	}
	GZ_CLOSE(reads_in);
	kseq_destroy(read);
	
	for (int i = 0; i < num_threads; ++i) {
		free_MappingThreadData(mdata[i]);
	}
	
	free_PackedDB(reference);
	destroy_lookup_table(lktbl);
	FCLOSE(out);
}
