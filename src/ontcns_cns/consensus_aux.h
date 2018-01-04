#ifndef CONSENSUS_AUX_H
#define CONSENSUS_AUX_H

#include "../common/gapped_candidate.h"
#include "../common/packed_db.h"
#include "../common/smart_assert.h"
#include "../common/output_stream.h"
#include "../gapped_align/ontcns_aligner.h"
#include "../tasc/cbcns.h"
#include "cns_options.h"

#include <set>

struct OverlapInfo
{
	int qid, qdir;
	int qoff, qend, qsize;
	int toff, tend, tsize;
	double ident_perc;
	size_t can_id;
	size_t align_size;
	size_t qalign_start;
	size_t talign_start;
};

struct OverlapInfoV
{
	std::vector<OverlapInfo> oiv;
	std::vector<char> align_string;
	
	void add(OntCnsAligner& aligner,
			 int qid_,
			 int qdir_,
			 int qsize_,
			 int tsize_,
			 size_t can_id_) {
		OverlapInfo oi;
		oi.qid = qid_;
		oi.qdir = qdir_;
		oi.qoff = aligner.query_start();
		oi.qend = aligner.query_end();
		oi.qsize = qsize_;
		oi.toff = aligner.target_start();
		oi.tend = aligner.target_end();
		oi.tsize = tsize_;
		oi.ident_perc = aligner.calc_ident_perc();
		oi.can_id = can_id_;
		oi.align_size = aligner.query_mapped_string().size();
		oi.qalign_start = align_string.size();
		oi.talign_start = oi.qalign_start + oi.align_size;
		oiv.push_back(oi);
		
		const std::string& qs = aligner.query_mapped_string();
		const std::string& ts = aligner.target_mapped_string();
		align_string.insert(align_string.end(), qs.begin(), qs.end());
		align_string.insert(align_string.end(), ts.begin(), ts.end());
	}
	
	void clear() {
		oiv.clear();
		align_string.clear();
	}
};

struct ConsensusData
{
	PackedDB*						reads;
	std::vector<GappedCandidate>*	candidates;
	std::size_t* 					next_can_id;
	pthread_mutex_t*				can_id_lock;
    CnsOptions*           			options;
    OutputStream               		out;
	OutputStream					raw_out;
	std::set<int>					corrected_read_ids;
	CntBaseConsensus				cbcns;
	OntCnsAligner					aligner;
	OntCnsAligner					aligner2;
	std::vector<int>				cov_stats;
	OverlapInfoV 					oiv_pb;
	OverlapInfoV					oiv_ont;

    ConsensusData(PackedDB* _reads,
				  std::vector<GappedCandidate>* _candidates,
				  std::size_t* _next_can_id,
				  pthread_mutex_t* _can_id_lock,
				  CnsOptions* _options,
				  std::ostream* _out,
				  pthread_mutex_t* _out_lock,
				  std::ostream* _raw_out,
				  pthread_mutex_t* _raw_out_lock) 
	:
	  reads(_reads),
	  candidates(_candidates),
	  next_can_id(_next_can_id),
	  can_id_lock(_can_id_lock),
	  options(_options),
	  out(*_out, _out_lock),
	  raw_out(*_raw_out, _raw_out_lock),
	  aligner(_options->error_rate, _options->align_method),
	  aligner2(_options->error_rate, _options->align_method)
	{ }
	
	bool next_cani_range(std::size_t& from, std::size_t& to) {
		bool r = 1;
		std::size_t i;
		int sid, sid1;
		pthread_mutex_lock(can_id_lock);
		i = *next_can_id;
		if (i >= candidates->size()) {
			r = 0;
		} else {
			from = i;
			sid = (*candidates)[i].sid;
			while (i != candidates->size()) {
				sid1 = (*candidates)[i].sid;
				if (sid != sid1) break;
				++i;
			}
			to = i;
			*next_can_id = i;
		}
		pthread_mutex_unlock(can_id_lock);
		return r;
	}
	
	~ConsensusData() {}
};

void
build_cns_thrd_data(PackedDB* _reads,
					std::vector<GappedCandidate>* _candidates,
					std::size_t* _next_can_id,
					pthread_mutex_t* _can_id_lock,
					CnsOptions* _options,
					std::ostream* _out,
					pthread_mutex_t* _out_lock,
					std::ostream* _raw_out,
					pthread_mutex_t* _raw_out_lock,
					ConsensusData** ppctd);

void
destroy_cns_thrd_data(ConsensusData** ppctd, const int num_threads);

#endif // CONSENSUS_AUX_H
