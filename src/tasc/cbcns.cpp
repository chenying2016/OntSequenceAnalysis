#include "cbcns.h"

#include "../common/cns_seq.h"

using namespace std;

void
normalise_gaps(const std::string& org_qaln, 
			   const std::string& org_taln,
			   std::string& qnorm,
			   std::string& tnorm)
{
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

void
CntBaseConsensus::add_one_align(const std::string& qaln, 
								const std::string& taln, 
								int toff,
								double weight,
							    bool need_normalise_gaps)
{
	if (need_normalise_gaps) {
		normalise_gaps(qaln, taln, nqaln, ntaln);
		get_cns_tags(nqaln, ntaln, toff, weight, tags);
	} else {
		get_cns_tags(qaln, taln, toff, weight, tags);
	}
}

bool
CntBaseConsensus::consensus(const int min_cov, 
							const int min_size, 
							const int template_id,
							const int template_size, 
							vector<pair<int, int> >& cns_intvs,
							OutputStream& out)
{
	build_backbone(tags, template_size, dci_alloc, li_alloc, backbone, coverage);
	int i = 0;
	CnsSeq cns_seq;
	while (i < template_size) {
		while (i < template_size && coverage[i] < min_cov) ++i;
		int j = i + 1;
		while (j < template_size && coverage[j] >= min_cov) ++j;
		if (j - i >= min_size * 0.85) {
			consensus_backbone_segment(backbone, i, j, coverage, cns_seq.seq);
			if (cns_seq.seq.size() >= min_size) {
				cns_seq.id = template_id;
				cns_seq.left = i;
				cns_seq.right = j;
				cns_seq.org_seq_size = template_size;
				out << cns_seq;
				cns_intvs.push_back(make_pair(i,j));
			}
		}
		i = j;
	}
	
	return cns_intvs.empty() ? 0 : 1;
}


struct CnsIntv
{
	int raw_from, raw_to;
	int cns_from, cns_to;
	
	CnsIntv(int a = 0, int b = 0, int c = 0, int d = 0): raw_from(a), raw_to(b), cns_from(c), cns_to(d) {}
};

bool
CntBaseConsensus::consensus(const int min_cov,
							const int min_size,
							const vector<char>& raw_read,
							string& cns_seq)
{
	cns_seq.clear();
	const int template_size = raw_read.size();
	build_backbone(tags, template_size, dci_alloc, li_alloc, backbone, coverage);
	int i = 0;
	int raw_from, raw_to, cns_from, cns_to;
	string cns_frag, cns_all;
	vector<CnsIntv> cns_intvs;
	while (i < template_size) {
		while (i < template_size && coverage[i] < min_cov) ++i;
		int j = i + 1;
		while (j < template_size && coverage[j] >= min_cov) ++j;
		if (j >= template_size) break;
		
		if (j - i >= min_size * 0.8) {
			consensus_backbone_segment(backbone, i, j, coverage, cns_frag, raw_from, raw_to);
			if (cns_frag.size() >= min_size) {
				cns_from = cns_all.size();
				cns_to = cns_from + cns_frag.size();
				cns_all.insert(cns_all.end(), cns_frag.begin(), cns_frag.end());
				cns_intvs.push_back(CnsIntv(raw_from, raw_to, cns_from, cns_to));
			}
		}
		i = j;
	}
	
	if (cns_intvs.empty()) {
		return 0;
	}
	
	int last_raw_to = 0;
	int last_cns_to = 0;
	for (size_t k = 0; k != cns_intvs.size(); ++k) {
		raw_from = cns_intvs[k].raw_from;
		raw_to = cns_intvs[k].raw_to;
		cns_from = cns_intvs[k].cns_from;
		cns_to = cns_intvs[k].cns_to;
		r_assert(cns_from == last_cns_to);
		if (raw_from > last_raw_to) {
            r_assert(raw_from < template_size);
			cns_seq.insert(cns_seq.end(), raw_read.begin() + last_raw_to, raw_read.begin() + raw_from);
		}
		cns_seq.insert(cns_seq.end(), cns_all.begin() + cns_from, cns_all.begin() + cns_to);
		last_raw_to = raw_to;
		last_cns_to = cns_to;
        r_assert(last_raw_to <= template_size);
	}
	
	if (last_raw_to != template_size) {
        r_assert(last_raw_to < template_size);
		cns_seq.insert(cns_seq.end(), raw_read.begin() + last_raw_to, raw_read.begin() + template_size);
	}
	
	for (size_t k = 0; k != cns_seq.size(); ++k) {
		int c = cns_seq[k];
		r_assert(c >= 0 && c < 4)(k)(c)(cns_seq.size());
		cns_seq[k] = "ACGT"[c];
	}
	
	return 1;
}
