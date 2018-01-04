#include "can_finder.h"
#include "mapping_aux.h"
#include "../word_finder/candidate_store.h"
#include "../common/timer.h"
#include "../common/pdb_aux.h"
#include "../word_finder/word_finder.h"

#include <fstream>

using namespace std;

void*
can_worker(void* arg)
{
	MappingData* md 				= (MappingData*)arg;
	OutputStream& out 				= md->out;
	PackedDB* reads 				= md->reads;
	PackedDB* reference 			= md->reference;
	int reads_start_id 				= md->reads_start_id;
	int reference_start_id 			= md->reference_start_id;
	LookupTable* lktbl 				= md->lktbl;
	InitWordParameters word_param 	= md->options->get_word_params();
	
	WordFinder word_finder(&word_param, lktbl, reference, reference_start_id);
	int sid, eid;
	vector<char> read_strands[2]; 
	read_strands[0].reserve(MAX_SEQ_SIZE);
	read_strands[1].reserve(MAX_SEQ_SIZE);
	vector<char> subject; 
	subject.reserve(MAX_SEQ_SIZE);
	CandidateStore can_store(word_param.num_candidates);
	
	while (md->get_next_chunk(sid, eid)) {
		for (int i = sid; i < eid; ++i) {
			int rsize = reads->seq_size(i);
			read_strands[FWD].resize(rsize);
			reads->get_sequence(i, true, read_strands[FWD].data());
			read_strands[REV].resize(rsize);
			reads->get_sequence(i, false, read_strands[REV].data());
			can_store.clear();
			word_finder.clear();
			word_finder.calc_candidates(read_strands[FWD].data(), i + reads_start_id, FWD, rsize, can_store, true); 
			word_finder.calc_candidates(read_strands[REV].data(), i + reads_start_id, REV, rsize, can_store, true);
			can_store.fix_sequence_ids(0, reference_start_id);
			can_store.output(out);
		}
	}
	
	return NULL;
}

idx
can_process_one_volume(MapOptions& options,
					  const idx reference_start_id,
					  const int svid,
					  const int evid,
					  const char* wrk_dir,
					  ostream& out)
{
	string vname;
	PackedDB reference;
	PackedDB reads;
	LookupTable lktbl;
	pthread_mutex_t out_lock;
	pthread_mutex_t chunk_lock;
	int chunk_id = 0;
	const int num_threads = options.num_threads;
	pthread_t tids[num_threads];
	MappingData* md[num_threads];
	idx reads_start_id = reference_start_id;
	
	make_volume_name(wrk_dir, svid, vname);
	reference.alloc(DEFAULT_VOLUME_SIZE + MAX_SEQ_SIZE);
	reference.load(vname.c_str());
	lktbl.build(&reference, options.kmer_size, options.kmer_cnt_cutoff, num_threads);
	reads.alloc(DEFAULT_VOLUME_SIZE + MAX_SEQ_SIZE);
	pthread_mutex_init(&out_lock, NULL);
	pthread_mutex_init(&chunk_lock, NULL);
	
	for (int i = 0; i < num_threads; ++i) {
		md[i] = new MappingData(i, &options, &reads, &reference, &lktbl, out, &out_lock, &chunk_id, &chunk_lock);
	}
	
	char process_msg[256];
	
	for (int vid = svid; vid < evid; ++vid) {
		sprintf(process_msg, "candidate detect volume %d vs volume %d.", vid, svid);
		DynamicTimer dt(process_msg);
		
		make_volume_name(wrk_dir, vid, vname);
		reads.load(vname.c_str());
		chunk_id = 0;
		for (int j = 0; j < num_threads; ++j) {
			md[j]->reset_reads_info(reads_start_id, reference_start_id);
			pthread_create(tids + j, NULL, can_worker, md[j]);
		}
		for (int j = 0; j < num_threads; ++j) {
			pthread_join(tids[j], NULL);
		}
		reads_start_id += reads.num_seqs();
	}
	
	for (int i = 0; i < num_threads; ++i) {
		delete md[i];
	}
	
	return reference.num_seqs();
}

void
candidate_main(MapOptions& options)
{
	const char* wrk_dir = options.wrk_dir;
	const int num_vols = load_num_volumes(wrk_dir);
	ostream* out = (options.output) ? (new ofstream(options.output)) : (&cout);
	idx reference_start_id = 0;
	for (int i = 0; i < num_vols; ++i) {
		reference_start_id += can_process_one_volume(options, reference_start_id, i, num_vols, wrk_dir, *out);
	}
	if (options.output) {
		delete out;
	}
}
