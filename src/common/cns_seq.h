#ifndef CNS_SEQ_H
#define CNS_SEQ_H

#include <string>
#include <iostream>

#include "output_stream.h"
#include "sequence.h"

struct CnsSeq
{
	int id;
	int left, right;
	int org_seq_size;
	std::string seq;
	
	template <class OStream>
	void write_to(OStream& os) const {
		std::string ontcns_hdr;
		make_ontcns_header(id, ontcns_hdr);
		os << '>'
		   << ontcns_hdr
		   << org_seq_size
		   << '/'
		   << left
		   << '_'
		   << right
		   << '_'
		   << seq.size()
		   << '\n'
		   << seq
		   << '\n';
	}
};

std::ostream& operator << (std::ostream& out, const CnsSeq& cns_seq);

OutputStream& operator << (OutputStream& out, const CnsSeq& cns_seq);

#endif // CNS_SEQ_H
