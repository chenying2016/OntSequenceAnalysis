#ifndef GAPPED_CANDIDATE_H
#define GAPPED_CANDIDATE_H

#include "aux_tools.h"
#include "defs.h"
#include "output_stream.h"

#include <iostream>

struct GappedCandidate
{
    int qid;
    int sid;
    int qdir;
    int sdir;
    int score;
	idx qbeg, qend, qsize;
	idx sbeg, send, ssize;
	idx qoff, soff;
	
	template <class OStream>
	void write_to(OStream& os) const {
		char buf[1024];
		char* p = buf;
		const char sep = '\t';
		int c;
		
		c = int2string(p, qid);   p += c; *p++ = sep;
		c = int2string(p, sid);   p += c; *p++ = sep;
		c = int2string(p, score); p += c; *p++ = sep;
		c = int2string(p, qdir);  p += c; *p++ = sep;
		c = int2string(p, qbeg);  p += c; *p++ = sep;
		c = int2string(p, qend);  p += c; *p++ = sep;
		c = int2string(p, qoff);  p += c; *p++ = sep;
		c = int2string(p, qsize); p += c; *p++ = sep;
		c = int2string(p, sdir);  p += c; *p++ = sep;
		c = int2string(p, sbeg);  p += c; *p++ = sep;
		c = int2string(p, send);  p += c; *p++ = sep;
		c = int2string(p, soff);  p += c; *p++ = sep;
		c = int2string(p, ssize); p += c; 
		
		*p++ = '\n';
		*p = '\0';
		os << buf;
	}
	
	template <class IStream>
	void read_from(IStream& is) {
		if (!(is >> qid)) return;
		  is >> sid
		   >> score
		   >> qdir
		   >> qbeg
		   >> qend
		   >> qoff
		   >> qsize
		   >> sdir
		   >> sbeg
		   >> send
		   >> soff
		   >> ssize;
	}
};

struct GappedCandidateScore_GT
{
    bool operator()(const GappedCandidate& a, const GappedCandidate& b) {
        return a.score > b.score;
    }
};

struct GappedCandidateSid_LT
{
    bool operator()(const GappedCandidate& a, const GappedCandidate& b) {
        return a.sid < b.sid;
    }
};

std::istream& operator>>(std::istream& in, GappedCandidate& can);

std::ostream& operator << (std::ostream& out, const GappedCandidate& m4);

OutputStream& operator << (OutputStream& out, const GappedCandidate& m4);

#endif // GAPPED_CANDIDATE_H
