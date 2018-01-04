#ifndef INDEL_DETECTION_SEQUENCE_H
#define INDEL_DETECTION_SEQUENCE_H

#include "buffer_line_reader.h"

#include <iostream>

class Sequence
{
public:
    typedef BufferLineReader LineReader;
    typedef BufferLineReader::OneDataLine str_t;

public:
    str_t& header() { return header_; }
    const str_t& header() const { return header_; }
    str_t& sequence() { return sequence_; }
    const str_t& sequence() const { return sequence_; }
    idx size() { return sequence_.size(); }
    idx size() const { return sequence_.size(); }
    void clear() { header_.clear(); sequence_.clear(); }

    idx read_one_seq(LineReader& line_reader);

private:
    str_t header_;
    str_t sequence_;
};

bool
is_ontcns_header(Sequence& seq);

void
make_ontcns_header(const int id, std::string& hdr);

int
extract_read_id_from_ontcns_hdr(const char* hdr);

void
check_and_rename_ontcns_header(Sequence& seq, const int id);

std::ostream& operator<<(std::ostream& out, const Sequence& seq);

#endif //INDEL_DETECTION_SEQUENCE_H
