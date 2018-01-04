#include "overlaps_provider.h"

using namespace std;

inline bool
check_mapping_ratio(GapAligner* aligner, int qsize, int tsize, double min_mapping_ratio)
{
    int qoff = aligner->query_start();
    int qend = aligner->query_end();
    int toff = aligner->target_start();
    int tend = aligner->target_end();

    return (qend - qoff >= qsize * min_mapping_ratio)
           ||
           (tend - toff >= tsize * min_mapping_ratio);
}

inline char
rc_char(const char c)
{
    switch (c) {
        case 'A':
            return 'T';
        case 'C':
            return 'C';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        default:
            cerr << "Error char " << c << "\n";
            exit(1);
    }
}

bool 
OverlapsProvider::extend_candidate(EGCInfo& cani, const double min_mapping_ratio, const int min_align_size)
{
    EGappedCandidate& can = *cani.p;
    if (can.status == EGappedCandidate::kInvalid) return 0;

    int qid, qdir, qext, qsize;
    int tid, tdir, text, tsize;
    vector<char>* qstr = NULL;
    vector<char>* tstr = NULL;

    if (cani.tid == can.qid) {
        qid     = can.sid;
        qdir    = can.sdir;
        qext    = can.sext;
        qsize   = can.ssize;
        qstr    = &raw_strands[2];
        tid     = can.qid;
        tdir    = can.qdir;
        text    = can.qext;
        tsize   = can.qsize;
        tstr    = &raw_strands[can.qdir];
    } else {
        qid     = can.qid;
        qdir    = can.qdir;
        qext    = can.qext;
        qsize   = can.qsize;
        qstr    = &raw_strands[can.qdir];
        tid     = can.sid;
        tdir    = can.sdir;
        text    = can.sext;
        tsize   = can.ssize;
        tstr    = &raw_strands[2];
    }

    if (can.status == EGappedCandidate::kAligned) {
        op.get_scripts(can, scripts);
		scripts_to_align_strings(*qstr, *tstr, can);
		can.status = EGappedCandidate::kInvalid;
		return 1;
    } else {
        bool r = aligner->go(qstr->data(), qext, qsize, tstr->data(), text, tsize, min_align_size);
        if (r && check_mapping_ratio(aligner, qsize, tsize, min_mapping_ratio)) {
            const string& qaln = aligner->query_mapped_string();
            const string& taln = aligner->target_mapped_string();
            build_edit_scripts(qaln, taln, scripts);
            op.add_scripts(scripts, can.es_start);
            if (can.es_start != -1) {
                can.es_size = (int)scripts.size();
                can.status = EGappedCandidate::kAligned;
                can.qoff = aligner->query_start();
                can.qend = aligner->query_end();
                can.soff = aligner->target_start();
                can.send = aligner->target_end();
            }
            true_qaln = qaln;
            true_taln = taln;
            return 1;
        }
    }

    return 0;
}

void 
OverlapsProvider::scripts_to_align_strings(std::vector<char>& query, std::vector<char>& target, EGappedCandidate& can)
{
    int qoff = can.qoff;
    int qend = can.qend;
    int toff = can.soff;
    int tend = can.send;
    const char* dt = "ACGT";
    const int OP_M = EditScriptWorker::kOpM;
    const int OP_I = EditScriptWorker::kOpI;
    const int OP_D = EditScriptWorker::kOpD;
    int type, num;
    true_qaln.clear();
    true_taln.clear();
    int qc, tc;

#define extract_qc qc = query[qoff]; r_assert(qc >= 0 && qc < 4)(qc); true_qaln += dt[qc];
#define extract_tc tc = target[qoff]; r_assert(tc >= 0 && tc < 4)(tc); true_taln += dt[tc];

    for (size_t k = 0; k != scripts.size(); ++k) {
        EditScriptWorker::unpack_edit_script(scripts[k], type, num);
        if (type != OP_M && type != OP_I && type != OP_D) {
            cerr << "k = " << k
                     << ", type = " << type
                     << "at " << can
                     << "\n";
            exit(1);
        }
        r_assert(num && num <= EditScriptWorker::kMaxNum)(num);
        switch (type) {
            case OP_M:
                for (int i = 0; i < num; ++i, ++qoff, ++toff) {
                    extract_qc;
                    extract_tc;
                    r_assert(qc == tc);
                }
                break;
            case OP_I:
                for (int i = 0; i < num; ++i, ++qoff) {
                    extract_qc;
                    true_taln += GAP_CHAR;
                }
                break;
            case OP_D:
                for (int i = 0; i < num; ++i, ++toff) {
                    true_qaln += GAP_CHAR;
                    extract_tc;
                }
                break;
            default:
                break;
        }
    }
    r_assert(qoff == qend)(qoff)(qend);
    r_assert(toff == tend)(toff)(tend);

#undef extract_qc
#undef extract_tc
}

void 
OverlapsProvider::adjust_overlap_info(EGCInfo&cani, EGappedCandidate& can)
{
    string* qaln_src;
    string* taln_src;
    bool reversed;
    const EGappedCandidate& csrc = *cani.p;
    if (cani.tid == csrc.qid) {
        qaln_src = &true_taln;
        taln_src = &true_qaln;
        if (csrc.qdir == FWD) {
            reversed = false;
            can.qoff = csrc.soff;
            can.qend = csrc.send;
            can.soff = csrc.qoff;
            can.send = csrc.qend;
        } else {
            reversed = true;
            can.qoff = csrc.ssize - csrc.send;
            can.qend = csrc.ssize - csrc.soff;
            can.soff = csrc.qsize - csrc.qend;
            can.send = csrc.qsize - csrc.qoff;
        }
    } else {
        qaln_src = &true_qaln;
        taln_src = &true_taln;
        if (csrc.sdir == FWD) {
            reversed = false;
            can.qoff = csrc.qoff;
            can.qend = csrc.qend;
            can.soff = csrc.soff;
            can.send = csrc.send;
        } else {
            reversed = true;
            can.qoff = csrc.qsize - csrc.qend;
            can.qend = csrc.qsize - csrc.qoff;
            can.soff = csrc.ssize - csrc.send;
            can.send = csrc.ssize - csrc.soff;
        }
    }

    if (reversed) {
        adjust_qaln.clear();
        adjust_taln.clear();
        size_t n = qaln_src->size();
        for (size_t i = n; i; --i) {
            char c = (*qaln_src)[i - 1];
			adjust_qaln += rc_char(c);
			c = (*taln_src)[i - 1];
			adjust_taln += rc_char(c);
        }
    } else {
        adjust_qaln = *qaln_src;
        adjust_taln = *taln_src;
    }
}

void
OverlapsProvider::normalise_gaps(const std::string& org_qaln, const std::string& org_taln)
{
	string& qnorm = normalised_qaln;
	string& tnorm = normalised_taln;
	qnorm.clear();
	tnorm.clear();
	
	r_assert(org_qaln.size() == org_taln.size());
	int qcnt = 0, tcnt = 0;
	for (size_t i = 0; i != org_qaln.size(); ++i) {
		if (org_qaln[i] != GAP_CHAR) ++qcnt;
		if (org_taln[i] != GAP_CHAR) ++tcnt;
	}
	
	for (size_t i = 0; i != org_qaln.size(); ++i) {
		const char qc = org_qaln[i];
		const char tc = org_taln[i];
		r_assert(qc != GAP_CHAR || tc != GAP_CHAR);
		if (qc != tc && qc != GAP_CHAR && tc != GAP_CHAR) {
			qnorm += GAP_CHAR;
			qnorm += qc;
			tnorm += tc;
			tnorm += GAP_CHAR;
		} else {
			qnorm += qc;
			tnorm += tc;
		}
	}
	
	size_t qlen = qnorm.size();
	size_t tlen = tnorm.size();
	for (size_t i = 0; i < qlen - 1; ++i) {
		// push target gaps
		if (tnorm[i] == GAP_CHAR) {
			size_t j = i;
			while (1) {
				const char c = tnorm[++j];
				if (c != GAP_CHAR || j > qlen - 1) {
					if (c == qnorm[i]) {
						tnorm[i] = c;
						tnorm[j] = GAP_CHAR;
					}
					break;
				}
			}
		}
		// push query gaps
		if (qnorm[i] == GAP_CHAR) {
			size_t j = i;
			while (1) {
				const char c = qnorm[++j];
				if (c != GAP_CHAR || j > tlen - 1) {
					if (c == tnorm[i]) {
						qnorm[i] = c;
						qnorm[j] = GAP_CHAR;
					}
					break;
				}
			}
		}
	}
	
	r_assert(qnorm.size() == tnorm.size());
	
	int qcnt2 = 0, tcnt2 = 0;
	for (size_t i = 0; i != qnorm.size(); ++i) {
		if (qnorm[i] != GAP_CHAR) ++qcnt2;
		if (tnorm[i] != GAP_CHAR) ++tcnt2;
	}
	r_assert(qcnt == qcnt2)(qcnt)(qcnt2);
	r_assert(tcnt == tcnt2);
}

bool 
OverlapsProvider::get_overlaps(EGCInfo& cani,
				  const int min_align_size,
				  const double min_mapping_ratio,
				  int& _qid,
				  int& qoff,
				  int& qend,
				  int& soff,
				  int& send,
				  const bool need_normalise_gaps)
{
	int qid, qdir, qext, qsize;
	int sid, sdir, sext, ssize;
	cani.extract_info(qid, qdir, qext, qsize, sid, sdir, sext, ssize, false);
	_qid = qid;
	r_assert(sid == cani.tid);
	r_assert(sdir == FWD);
	r_assert(ssize == raw_strands[FWD].size());
	vector<char>& qstr = raw_strands[2];
	vector<char>& tstr = raw_strands[FWD];
	qstr.resize(qsize);
	reads.get_sequence(qid, qdir == FWD, qstr.data());
	bool r = aligner->go(qstr.data(), qext, qsize, tstr.data(), sext, ssize, min_align_size);
	if (r) r = check_mapping_ratio(aligner, qsize, ssize, min_mapping_ratio);
	if (!r) return 0;
	if (need_normalise_gaps) {
		normalise_gaps(aligner->query_mapped_string(), aligner->target_mapped_string());
	} else {
		normalised_qaln = aligner->query_mapped_string();
		normalised_taln = aligner->target_mapped_string();
	}
	qoff = aligner->query_start();
	qend = aligner->query_end();
	soff = aligner->target_start();
	send = aligner->target_end();
	
	///*
	int a = qoff, b = soff;
	const string& aaln = aligner->query_mapped_string();
	const string& baln = aligner->target_mapped_string();
	r_assert(aaln[0] == baln[0])(aaln[0])(baln[0]);
	r_assert(aaln[ aaln.size() - 1 ] == baln[ baln.size() - 1]);
	for (size_t i = 0; i != aaln.size(); ++i) {
		char qc = aaln[i];
		if (qc != '-') {
			char qc1 = qstr[a];
			qc1 = "ACGT"[qc1];
			r_assert(qc == qc1)(a)(i)(qc)(qc1);
			++a;
		}
		char sc = baln[i];
		if (sc != '-') {
			char sc1 = tstr[b];
			sc1 = "ACGT"[sc1];
			r_assert(sc == sc1);
			++b;
		}
	}
	
	return 1;
}

int 
OverlapsProvider::global_read_id(const int local_id) const
{
	static int cnt = 0;
	const char* hdr = reads.seq_header(local_id);
	int gid = atoi(hdr);
	if (!cnt) {
		cout << hdr << ' ' << gid << endl;
		cnt = 1;
	}
	return gid;
}
