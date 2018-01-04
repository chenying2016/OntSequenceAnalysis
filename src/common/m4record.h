#ifndef M4RECORD_H
#define M4RECORD_H

#include "aux_tools.h"
#include "defs.h"
#include "smart_assert.h"
#include "output_stream.h"

#include <cstdio>
#include <iostream>

/* There are two set of mapped ranges.
 * Let [dir, begin, end, size] be the mapped range of sequence S and S' be the reversed complement of S.
 * 1) The normal range. The offsets are relateve to the forward strand always.
 *    If dir == FWD, forward strand, then S[begin, end) is mapped.
 *    If dir == REV, reversed strand, then S'[size - end, size - begin) is mapped.
 * 2) The true range. The offsets are set according to the strand.
 *    If dir == FWD, forward strand, then S[begin, end) is mapped.
 *    If dir == REV, reversed strand, then S'[begin, end) is mapped.
 */

#define invalid_ext (-1)

struct M4Record
{
    idx qid;
    int qdir;
    idx qoff;
    idx qend;
    idx qext;
    idx qsize;
    idx sid;
    int sdir;
    idx soff;
    idx send;
    idx sext;
    idx ssize;
    double ident_perc;
    int vscore;
	
	template <class OStream>
	void write_to(OStream& os) const {
		char buf[1024];
		char* p = buf;
		const char sep = '\t';
		int c;
		
		c = int2string(p, qid);      p += c; *p++ = sep;
		c = int2string(p, sid);      p += c; *p++ = sep;
		c = snprintf(p, 32, "%.2lf", ident_perc); 
									 p += c; *p++ = sep;
		c = int2string(p, vscore);   p += c; *p++ = sep;
		c = int2string(p, qdir);     p += c; *p++ = sep;
		c = int2string(p, qoff);     p += c; *p++ = sep;
		c = int2string(p, qend);     p += c; *p++ = sep;
		c = int2string(p, qsize);    p += c; *p++ = sep;
		c = int2string(p, sdir);     p += c; *p++ = sep;
		c = int2string(p, soff);     p += c; *p++ = sep;
		c = int2string(p, send);     p += c; *p++ = sep;
		c = int2string(p, ssize);    p += c;
		
		if (qext != invalid_ext) {
			*p++ = sep;
			c = int2string(p, qext); p += c; *p++ = sep;
			c = int2string(p, sext); p += c; 
		}
		
		*p++ = '\n';
		*p = '\0';
		os << buf;
	}
	
	template <class IStream>
	void read_from(IStream& is) {
		if (!(is >> qid)) return;
		is >> sid
		   >> ident_perc
		   >> vscore
		   >> qdir
		   >> qoff
		   >> qend
		   >> qsize
		   >> sdir
		   >> soff
		   >> send
		   >> ssize;
		if (is >> qext) is >> sext;
	}
};

inline idx m4qid(const M4Record& m) { return m.qid; }
inline idx m4sid(const M4Record& m) { return m.sid; }
inline int m4qdir(const M4Record& m) { return m.qdir; }
inline int m4sdir(const M4Record& m) { return m.sdir; }
inline idx m4qs(const M4Record& m) { return m.qsize; }
inline idx m4ss(const M4Record& m) { return m.ssize; }
inline double m4ident(const M4Record& m) { return m.ident_perc; }
inline int m4vscore(const M4Record& m) { return m.vscore; }
inline int m4size(const M4Record& m) { return m.qend - m.qoff; }

inline idx m4nqb(const M4Record& m) { return m.qoff; }
inline idx m4nqe(const M4Record& m) { return m.qend; }
inline idx m4nqg(const M4Record& m) { return m.qext; }
inline idx m4tqb(const M4Record& m) { return (m.qdir == FWD) ? (m.qoff) : (m.qsize - m.qend); }
inline idx m4tqe(const M4Record& m) { return (m.qdir == FWD) ? (m.qend) : (m.qsize - m.qoff); }
inline idx m4tqg(const M4Record& m) { return (m.qdir == FWD) ? (m.qext) : (m.qsize - 1 - m.qext); }

inline idx m4nsb(const M4Record& m) { return m.soff; }
inline idx m4nse(const M4Record& m) { return m.send; }
inline idx m4nsg(const M4Record& m) { return m.sext; }
inline idx m4tsb(const M4Record& m) { return (m.sdir == FWD) ? (m.soff) : (m.ssize - m.send); }
inline idx m4tse(const M4Record& m) { return (m.sdir == FWD) ? (m.send) : (m.ssize - m.soff); }
inline idx m4tsg(const M4Record& m) { return (m.sdir == FWD) ? (m.sext) : (m.ssize - 1 - m.sext); }

inline void set_m4qid(const idx qid, M4Record& m) { m.qid = qid; }
inline void set_m4sid(const idx sid, M4Record& m) { m.sid = sid; }

inline void
fill_m4_mapped_range(idx id, int dir, idx off, idx end, idx ext, idx size,
        bool fill_query_range, bool from_normal_range, M4Record& m)
{
    idx *pid, *poff, *pend, *pext, *psize;
    int *pdir;
    if (fill_query_range) {
        pid     = &m.qid;
        pdir    = &m.qdir;
        poff    = &m.qoff;
        pend    = &m.qend;
        pext    = &m.qext;
        psize   = &m.qsize;
    } else {
        pid     = &m.sid;
        pdir    = &m.sdir;
        poff    = &m.soff;
        pend    = &m.send;
        pext    = &m.sext;
        psize   = &m.ssize;
    }

    *pid    = id;
    *pdir   = dir;
    *poff   = off;
    *pend   = end;
    *pext   = ext;
    *psize  = size;

    if ((!from_normal_range) && (dir == REV)) {
        *poff = size - end;
        *pend = size - off;
        *pext = size - 1 - ext;
    }
	
	r_assert(*poff >= 0);
}

inline void
fill_m4record(idx qid, int qdir, idx qoff, idx qend, idx qext, idx qsize,
              idx sid, int sdir, idx soff, idx send, idx sext, idx ssize,
              double ident_perc, int vscore, bool from_normal_range,
              M4Record& m)
{
    fill_m4_mapped_range(qid, qdir, qoff, qend, qext, qsize, true, from_normal_range, m);
    fill_m4_mapped_range(sid, sdir, soff, send, sext, ssize, false, from_normal_range, m);
    m.ident_perc = ident_perc;
    m.vscore = vscore;
}

struct M4RecordSize_GT
{
	bool operator()(const M4Record& a, const M4Record& b) const {
		return m4size(a) > m4size(b);
	}
};

inline bool
is_full_map_ratio(const M4Record& m)
{
	return (m4nqe(m) - m4nqb(m) >= m4qs(m) * 0.9)
		   ||
		   (m4nse(m) - m4nsb(m) >= m4ss(m) * 0.9);
}

inline bool
is_full_map_range(const M4Record& m, const int unmapped_end_size = 500)
{
	return (m.qoff < unmapped_end_size && m.qsize - m.qend < unmapped_end_size)
		   ||
		   (m.soff < unmapped_end_size && m.ssize - m.send < unmapped_end_size);
}

inline bool
is_left_child(const M4Record& host, const M4Record& slave)
{
    if (m4qdir(host) != m4qdir(slave)) return 0;
    if (m4sid(host) != m4sid(slave)) return 0;
    
    const idx hqb = m4tqb(host);
    const idx hsb = m4tsb(host);
    const idx sqe = m4tqe(slave);
    const idx sse = m4tse(slave);
	
	const int close_dist = 200;

    /// read close
    if (abs(hqb - sqe) <= close_dist) return 1;

    /// reference close
    if (abs(hsb - sse) <= close_dist) return 1;

    return 0;
}

inline bool
is_right_child(const M4Record& host, const M4Record& slave)
{
    if (m4qdir(host) != m4qdir(slave)) return 0;
    if (m4sid(host) != m4sid(slave)) return 0;

    const idx hqe = m4tqe(host);
    const idx hse = m4tse(host);
    const idx sqb = m4tqb(slave);
    const idx ssb = m4tsb(slave);
	
	const int close_dist = 200;

    /// read close
    if (abs(hqe - sqb) <= close_dist) return 1;

    /// reference close
    if (abs(hse - ssb) <= close_dist) return 1;

    return 0;
}

std::istream& operator >> (std::istream& in, M4Record& m);

std::ostream& operator << (std::ostream& out, const M4Record& m);

OutputStream& operator << (OutputStream& out, const M4Record& m4);

#endif // M4RECORD_H
