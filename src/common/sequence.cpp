#include "sequence.h"

#include <sstream>
#include <string>

#include "mc_log.h"

using namespace std;

static const string kOntCnsHeader = "ontcns/";

std::ostream& operator<<(std::ostream& out, const Sequence& seq)
{
    out << ">";
    const Sequence::str_t& header = seq.header();
    idx n = header.size();
    for (idx i = 0; i < n; ++i) out << header[i];
    out << "\n";

    const Sequence::str_t& s = seq.sequence();
    n = s.size();
    for (idx i = 0; i < n; ++i) out << s[i];
    out << "\n";

    return out;
}

bool
is_ontcns_header(Sequence& seq)
{
	Sequence::str_t& hdr = seq.header();
	if (hdr.size() < kOntCnsHeader.size()) return 0;

	for (size_t i = 0; i != kOntCnsHeader.size(); ++i) {
		if (hdr[i] != kOntCnsHeader[i]) return 0;
	}
	return 1;
}

void
make_ontcns_header(const int id, std::string& hdr)
{
	ostringstream os;
	os << id << '/';
	hdr = kOntCnsHeader;
	hdr += os.str();
}

void
check_and_rename_ontcns_header(Sequence& seq, const int id)
{
	if (is_ontcns_header(seq)) return;
	string new_hdr;
	make_ontcns_header(id, new_hdr);
	Sequence::str_t& hdr = seq.header();
	hdr.assign(new_hdr.begin(), new_hdr.end());
}

int
extract_read_id_from_ontcns_hdr(const char* hdr)
{
	const char* p = hdr;
	while (*p != '\0' && *p != '/') ++p;
	if (*p == '\0') {
		mc_error << "invalid ontcns header: " << hdr << eolog;
	}
	++p;
	return atoi(p);
}
