#ifndef ONTCNS_ALIGNER_H
#define ONTCNS_ALIGNER_H

#include <string>
#include <vector>

#include "diff_ex.h"
#include "edlib_ex.h"
#include "xdrop_sw.h"

#include "../common/mc_log.h"

#define ALIGN_METHOD_EDLIB 0
#define ALIGN_METHOD_DIFF  1

class OntCnsAligner
{
public:
    OntCnsAligner(const double e, int align_method = ALIGN_METHOD_EDLIB) {
		edlib_aligner = NULL;
		diff_aligner = NULL;
		if (align_method == ALIGN_METHOD_EDLIB) {
			edlib_aligner = new EdLibAligner;
		} else if (align_method == ALIGN_METHOD_DIFF) {
			diff_aligner = new DiffAligner;
		} else {
			mc_error << "Invalid Method: " << align_method << eolog;
		}
        error = e;
        qabuf = new char[100000];
        tabuf = new char[100000];
    }

    ~OntCnsAligner() {
		if (edlib_aligner) {
			delete edlib_aligner;
		}
		if (diff_aligner) {
			delete diff_aligner;
		}
        delete[] qabuf;
        delete[] tabuf;
    }

    double calc_ident_perc() const {
        return ident_perc;
    }

    int query_start() const {
        return qoff;
    }

    int query_end() const {
        return qend;
    }

    int target_start() const {
        return toff;
    }

    int target_end() const {
        return tend;
    }
	
	void set_ids(int _qid, int _tid) {
		qid = _qid;
		tid = _tid;
	}

    const std::string& query_mapped_string() const {
        return query_align;
    }

    const std::string& target_mapped_string() const {
        return target_align;
    }

    bool go(const char* query, const int query_start, const int query_size,
            const char* target, const int target_start, const int target_size,
            const int min_align_size);

private:
    void align_ex(const char* query,
            const int query_size,
            const char* target,
            const int target_size,
            const bool forward_extend,
            std::string& qaln,
            std::string& taln);

private:
    EdLibAligner*       edlib_aligner;
    //XdropSwAligner      xdrop_aligner;
	DiffAligner*		diff_aligner;
    int                 qoff;
    int                 qend;
    int                 toff; 
    int                 tend;
    double              ident_perc;
    std::string         query_align;
    std::string         target_align;
    std::string         fqaln;
    std::string         rqaln;
    std::string         ftaln;
    std::string         rtaln;
    std::vector<char>   qfrag;
    std::vector<char>   tfrag;
    char*               qabuf;
    char*               tabuf;
    double              error;
	
	int					qid;
	int					tid;
};

#endif // ONTCNS_ALIGNER_H
