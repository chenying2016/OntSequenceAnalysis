#ifndef CANDIDATE_INFO_H
#define CANDIDATE_INFO_H

#include "../word_finder/word_finder_aux.h"

struct EGappedCandidate
{
    int qdir, qid, qext, qsize;
    int sdir, sid, sext, ssize;
    int score;
    int status;
    int qoff, qend;
    int soff, send;
    int es_start, es_size;

    void init(const GappedCandidate& c) {
        qdir        = c.qdir;
        qid         = c.qid;
        qext        = c.qoff;
        qsize       = c.qsize;
        sdir        = FWD;
        sid         = c.sid;
        sext        = c.soff;
        ssize       = c.ssize;
        score       = c.score;
        status      = kValid;
        es_start    = -1;
    }

    static const int kValid = 0;
    static const int kInvalid = 1;
    static const int kAligned = 2;
};

std::ostream& operator << (std::ostream& out, const EGappedCandidate& can);

struct EGCInfo
{
    int tid;
	int idx_from, idx_to;
    EGappedCandidate* p;

    int score() const {
        return p->score;
    }
	
	inline int template_size() const {
		return (tid == p->qid) ? (p->qsize) : (p->ssize);
	}
	
	inline int template_id() const {
		return tid;
	}
	
	inline int non_template_id() const {
		return (tid == p->qid) ? (p->sid) : (p->qid);
	}
	
	inline int non_template_size() const {
		return (tid == p->qid) ? (p->ssize) : (p->qsize);
	}
	
	void extract_info(int &qid,
					  int& qdir,
					  int& qext,
					  int& qsize,
					  int& sid,
					  int& sdir,
					  int& sext,
					  int& ssize,
					  bool true_offset) {
		if (tid == p->qid) {
			qid = p->sid;
			qsize = p->ssize;
			sid = p->qid;
			ssize = p->qsize;
			if ((!true_offset) && (p->qdir == REV)) {
				qdir = REVERSE_STRAND(p->sdir);
				qext = p->ssize - 1 - p->sext;
				sdir = REVERSE_STRAND(p->qdir);
				sext = p->qsize - 1 - p->qext;
			} else {
				qdir = p->sdir;
				qext = p->sext;
				sdir = p->qdir;
				sext = p->qext;
			}
		} else {
			qid = p->qid;
			qsize = p->qsize;
			sid = p->sid;
			ssize = p->ssize;
			if ((!true_offset) && (p->sdir == REV)) {
				qdir = REVERSE_STRAND(p->qdir);
				qext = p->qsize - 1 - p->qext;
				sdir = REVERSE_STRAND(p->sdir);
				sext = p->ssize - 1 - p->sext;
			} else {
				qdir = p->qdir;
				qext = p->qext;
				sdir = p->sdir;
				sext = p->sext;
			}
		}
	}
	
	void get_non_template_info(bool true_offset, int& id, int& ext, int& size) const {
		if (tid == p->qid) {
			id   = p->sid;
			ext  = p->sext;
			size = p->ssize;
			if ((!true_offset) && (p->sdir == REV)) ext = size - 1 - ext;
		} else {
			id 	 = p->qid;
			ext	 = p->qext;
			size = p->qsize;
			if ((!true_offset) && (p->qdir == REV)) ext = size - 1 - ext;
		}
	}
};

struct CmpEGCInfo_Score
{
    bool operator()(const EGCInfo& a, const EGCInfo& b) const {
        int ascore = a.score();
        int bscore = b.score();
        if (ascore != bscore) return ascore > bscore;
        int aid, aext, asize;
        a.get_non_template_info(false, aid, aext, asize);
        int bid, bext, bsize;
        b.get_non_template_info(false, bid, bext, bsize);
        if (aid != bid) return aid < bid;
        return aext < bext;
    }
};

struct CmpEGCInfo_Tid
{
	bool operator()(const EGCInfo& a, const EGCInfo& b) const {
		return a.tid < b.tid;
	}
};

#endif // CANDIDATE_INFO_H
