#include "ontcns_aligner.h"

#include "../common/smart_assert.h"

using namespace std;

static const int kMatCnt = 8;
static const int kBlockSize = 512;

bool
retrieve_next_sequence_block(const char* query,
        int qidx,
        const int qsize,
        const char* target,
        int tidx,
        const int tsize,
        const int desired_block_size,
        const bool forward_extend,
        vector<char>& qfrag,
        vector<char>& tfrag)
{
    bool last_block = 0;
    int qleft = qsize - qidx;
    int tleft = tsize - tidx;
    int qblk;
    int tblk;
    if (qleft < desired_block_size + 100 || tleft < desired_block_size + 100) {
        qblk = min(qleft, static_cast<int>(tleft * 1.3));
        tblk = min(tleft, static_cast<int>(qleft * 1.3));
        last_block = 1;
    } else {
        qblk = desired_block_size;
        tblk = desired_block_size;
        last_block = 0;
    }

    if (forward_extend) {
        const char* Q = query + qidx;
        qfrag.clear();
        for (int i = 0; i < qblk; ++i) qfrag.push_back( Q[i] );

        const char* R = target + tidx;
        tfrag.clear();
        for (int i = 0; i < tblk; ++i) tfrag.push_back( R[i] );
    } else {
        const char* Q = query - qidx;
        qfrag.clear();
        for (int i = 0; i < qblk; ++i) qfrag.push_back( Q[-i] );

        const char* R = target - tidx;
        tfrag.clear();
        for (int i = 0; i < tblk; ++i) tfrag.push_back( R[-i] );
    }

    return last_block;
}
    
void 
OntCnsAligner::align_ex(const char* query,
            const int query_size,
            const char* target,
            const int target_size,
            const bool forward_extend,
            string& qaln,
            string& taln)
{
	if (0)
	cout << "query_size = " << query_size
		 << ", target_size = " << target_size
		 << "\n";
    int qidx = 0, tidx = 0;
    qaln.clear();
    taln.clear();
    while (1) {
        int qfae, tfae, qfrag_size, tfrag_size;
        bool last_block = retrieve_next_sequence_block(query,
                qidx,
                query_size,
                target,
                tidx,
                target_size,
                kBlockSize,
                forward_extend,
                qfrag,
                tfrag);
        qfrag_size = qfrag.size();
        tfrag_size = tfrag.size();
        if (qfrag_size == 0 || tfrag_size == 0) break;

		if (edlib_aligner) {
			edlib_aligner->align(qfrag.data(),
                qfrag_size,
                tfrag.data(),
                tfrag_size,
                error,
                qabuf,
                tabuf,
                qfae,
                tfae);
		} else if (diff_aligner) {
			diff_aligner->align(qfrag.data(),
								 qfrag_size,
								 tfrag.data(),
								 tfrag_size,
								 error,
								 qabuf,
								 tabuf,
								 qfae,
								 tfae);
		}

        bool done = last_block;
        int acnt = 0, qcnt = 0, tcnt = 0;
        if (qfrag_size - qfae > 30 && tfrag_size - tfae > 30) done = 1;
        int align_size = strlen(qabuf);
        int k = align_size - 1, m = 0;
        while (k >= 0) {
            const char qc = qabuf[k];
            const char tc = tabuf[k];
            if (qc != GAP_CHAR) ++qcnt;
            if (tc != GAP_CHAR) ++tcnt;
            if (qc == tc) {
                ++m;
            } else {
                m = 0;
            }
            ++acnt;
            if (m == kMatCnt) break;
            --k;
        }

        if (m != kMatCnt || k < 1) {
            align_size = 0;
            for (int i = 0; i < qfrag_size && i < tfrag_size; ++i) {
                const char qc = qfrag[i];
                const char tc = tfrag[i];
                if (qc != tc) break;
                qabuf[align_size] = DNADecoder::decode(qc);
                tabuf[align_size] = DNADecoder::decode(tc);
                ++align_size;
            }
            done = 1;
        } else {
            align_size -= acnt;
            qidx += (qfae - qcnt);
            tidx += (tfae - tcnt);
            if (done) align_size += kMatCnt;
        }
		
		if (0)
		cout << "qidx = " << qidx
			 << ", tidx = " << tidx
			 << ", qfrag_size = " << qfrag_size
			 << ", qfae = " << qfae
			 << ", tfrag_size = " << tfrag_size
			 << ", tfae = " << tfae
			 << ", k = " << k
			 << ", m = " << m
			 << ", done = " << done
			 << "\n";

        qabuf[align_size] = '\0';
        tabuf[align_size] = '\0';
        qaln.insert(qaln.end(), qabuf, qabuf + align_size);
        taln.insert(taln.end(), tabuf, tabuf + align_size);
        if (done) break;
    }
}

bool
OntCnsAligner::go(const char* query,
        const int query_start,
        const int query_size,
        const char* target,
        const int target_start,
        const int target_size,
        const int min_align_size)
{
	//if (!(qid == 4277 && tid == 1)) return 1;
    int X = query_start;
    int XS = query_size;
    int Y = target_start;
    int YS = target_size;
    int fqcnt = 0, rqcnt = 0, ftcnt = 0, rtcnt = 0;
    query_align.clear();
    target_align.clear();

    align_ex(query + X - 1,
            X,
            target + Y - 1,
            Y,
            0,
            rqaln,
            rtaln);

    /// trim non-match prolog
    {
        int qcnt = 0, tcnt = 0, acnt = 0, m = 0;
        r_assert(rqaln.size() == rtaln.size());
        for (size_t i = 0; i != rqaln.size(); ++i) {
            const char qc = rqaln[i];
            const char tc = rtaln[i];
            if (qc != GAP_CHAR) ++qcnt;
            if (tc != GAP_CHAR) ++tcnt;
            if (qc == tc) {
                ++m;
            } else {
                m = 0;
            }
            ++acnt;
            if (m == kMatCnt) break;
        }
        if (m == kMatCnt) {
            X -= qcnt;
            Y -= tcnt;
            int n = rqaln.size();
            for (int i = n; i > acnt; --i) {
                const char qc = rqaln[i - 1];
                query_align.push_back(qc);
                if (qc != GAP_CHAR) ++rqcnt;
                const char tc = rtaln[i - 1];
                target_align.push_back(tc);
                if (tc != GAP_CHAR) ++rtcnt;
            }
        }
    }

    align_ex(query + X,
             XS - X,
             target + Y,
             YS - Y,
             1,
             fqaln,
             ftaln);

    /// trime non-match prolog
    if (query_align.empty()) {
        int qcnt = 0, tcnt = 0, acnt = 0, m = 0;
        r_assert(fqaln.size() == ftaln.size());
        for (size_t i = 0; i != fqaln.size(); ++i) {
            const char qc = fqaln[i];
            const char tc = ftaln[i];
            if (qc != GAP_CHAR) ++qcnt;
            if (tc != GAP_CHAR) ++tcnt;
            if (qc == tc) {
                ++m;
            } else {
                m = 0;
            }
            ++acnt;
            if (m == kMatCnt) break;
        }
        if (m == kMatCnt) {
            acnt -= kMatCnt;
            qcnt -= kMatCnt;
            tcnt -= kMatCnt;
            X += qcnt;
            Y += tcnt;
            int n = fqaln.size();
            for (int i = acnt; i < n; ++i) {
                const char qc = fqaln[i];
                if (qc != GAP_CHAR) ++fqcnt;
                query_align.push_back(qc);
                const char tc = ftaln[i];
                if (tc != GAP_CHAR) ++ftcnt;
                target_align.push_back(tc);
            }
        }
    } else {
        for (size_t i = 0; i != fqaln.size(); ++i) {
            const char qc = fqaln[i];
            if (qc != GAP_CHAR) ++fqcnt;
            query_align.push_back(qc);
            const char tc = ftaln[i];
            if (tc != GAP_CHAR) ++ftcnt;
            target_align.push_back(tc);
        }
    }

    qoff = X - rqcnt;
    qend = X + fqcnt;
    toff = Y - rtcnt;
    tend = Y + ftcnt;

    if (1) {
        int x = qoff, y = toff;
        for (size_t i = 0; i != query_align.size(); ++i) {
            const char qc = query_align[i];
            if (qc != GAP_CHAR) {
                const char qc1 = DNADecoder::decode(query[x]);
                r_assert(qc == qc1)(i)(x)(y)(qc)(qc1)(qid)(tid);
                ++x;
            }
            const char tc = target_align[i];
            if (tc != GAP_CHAR) {
                const char tc1 = DNADecoder::decode(target[y]);
                r_assert(tc == tc1)(i)(x)(y)(tc)(tc1)(qid)(tid);
                ++y;
            }
        }
    }

    r_assert(qoff >= 0)(qoff)(rqcnt)(X)(query_size)(query_start)(query_size)(target_start)(target_size);
    r_assert(qend <= query_size)(qend)(fqcnt)(X)(query_size);
    r_assert(toff >= 0)(toff)(rtcnt)(Y)(target_size);
    r_assert(tend <= target_size)(tend)(ftcnt)(Y)(target_size);

    int m = query_align.size(), n = 0;
    for (size_t i = 0; i != query_align.size(); ++i) {
        if (query_align[i] == target_align[i]) ++n;
    }
    if (m) {
        ident_perc = 100.0 * n / m;
    } else {
        ident_perc = 0.0;
    }

    return qend - qoff >= min_align_size || tend - toff >= min_align_size;
}
