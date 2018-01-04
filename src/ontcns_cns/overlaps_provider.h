#ifndef OVERLAPS_PROVIDER_H
#define OVERLAPS_PROVIDER_H

#include "overlap_pool.h"
#include "candidate_info.h"
#include "../gapped_align/ontcns_aligner.h"

class OverlapsProvider
{
public:
    OverlapsProvider(PackedDB& _reads):
        reads(_reads),
		aligner(new GapAligner(0.5)) {
            std::size_t s = MAX_SEQ_SIZE;
            scripts.reserve(s);
            raw_strands[0].reserve(s);
            raw_strands[1].reserve(s);
            raw_strands[2].reserve(s);
            true_qaln.reserve(s);
            true_taln.reserve(s);
            adjust_qaln.reserve(s);
            adjust_taln.reserve(s);
        }

    ~OverlapsProvider() {
        delete aligner;
    }

    void prefetch_target_sequences(const int tid) {
        std::size_t tsize = (std::size_t)reads.seq_size(tid);
        raw_strands[FWD].resize(tsize);
        raw_strands[REV].resize(tsize);
        reads.get_sequence(tid, true, raw_strands[FWD].data());
        reads.get_sequence(tid, true, raw_strands[REV].data());
    }
	
	const std::vector<char>& target_sequence() {
		return raw_strands[FWD];
	}
	
	const char* seq_header(const int id) const {
		return reads.seq_header(id);
	}
	
	bool get_overlaps(EGCInfo& cani,
					  const int min_align_size,
					  const double min_mapping_range,
					  int& qid,
					  int& qoff,
					  int& qend,
					  int& soff,
					  int& send,
					  const bool need_normalise_gaps);
	
	const std::string& nqaln() const {
		return normalised_qaln;
	}
	
	const std::string& ntaln() const {
		return normalised_taln;
	}
	
	int global_read_id(const int local_id) const;
	
	void set_min_ident(const double ident) {
		//aligner->set_min_ident(ident);
	}

private:
    bool extend_candidate(EGCInfo& cani, const double min_mapping_ratio, const int min_align_size);

    void scripts_to_align_strings(std::vector<char>& query, std::vector<char>& target, EGappedCandidate& can);

    void adjust_overlap_info(EGCInfo&cani, EGappedCandidate& can);
	
	void normalise_gaps(const std::string& org_qaln, const std::string& org_taln);

private:
    OverlapPool         op;
    std::vector<u8>     scripts;
    std::vector<char>   raw_strands[3];
    std::string         true_qaln;
    std::string         true_taln;
    std::string         adjust_qaln;
    std::string         adjust_taln;
	std::string			normalised_qaln;
	std::string			normalised_taln;
    PackedDB&           reads;
    GapAligner*         aligner;
};

#endif // OVERLAPS_PROVIDER_H
