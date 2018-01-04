#include "cns_seq.h"

std::ostream& operator << (std::ostream& out, const CnsSeq& cns_seq)
{
	cns_seq.write_to(out);
	return out;
}

OutputStream& operator << (OutputStream& out, const CnsSeq& cns_seq)
{
	MCOstream& mcos = out.mcostream();
	cns_seq.write_to(mcos);
	out.absorp_os_data();
	return out;
}
