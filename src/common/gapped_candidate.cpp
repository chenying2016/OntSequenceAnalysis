#include "gapped_candidate.h"

#include <sstream>

using namespace std;

std::istream&  operator>>(std::istream& in, GappedCandidate& can)
{
	string line;
	if (getline(in, line)) {
		istringstream ins(line);
		ins >> can.qid
			>> can.sid
			>> can.score
			>> can.qdir
			>> can.qbeg
			>> can.qend
			>> can.qoff
			>> can.qsize
			>> can.sdir
			>> can.sbeg
			>> can.send
			>> can.soff
			>> can.ssize;
	}
	return in;
}

std::ostream& operator << (std::ostream& out, const GappedCandidate& m4)
{
	m4.write_to(out);
	return out;
}

OutputStream& operator << (OutputStream& out, const GappedCandidate& m4)
{
	MCOstream& mcos = out.mcostream();
	m4.write_to(mcos);
	out.absorp_os_data();
	return out;
}
