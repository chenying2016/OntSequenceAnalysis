#include "m4record.h"

#include <sstream>
#include <string>

using namespace std;

istream& operator >> (istream& in, M4Record& m)
{
    string line;
    if (!getline(in, line)) return in;
    istringstream ins(line);
    ins >> m.qid
        >> m.sid
        >> m.ident_perc
        >> m.vscore
        >> m.qdir
        >> m.qoff
        >> m.qend
        >> m.qsize
        >> m.sdir
        >> m.soff
        >> m.send
        >> m.ssize;
    if (ins >> m.qext) ins >> m.sext;

    return in;
}

std::ostream& operator << (std::ostream& out, const M4Record& m)
{
	m.write_to(out);
	return out;
}

OutputStream& operator << (OutputStream& out, const M4Record& m4)
{
	MCOstream& mcos = out.mcostream();
	m4.write_to(mcos);
	out.absorp_os_data();
	return out;
}
