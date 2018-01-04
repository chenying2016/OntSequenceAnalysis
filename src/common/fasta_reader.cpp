//
// Created by sysu on 2/26/17.
//

#include "fasta_reader.h"
#include "aux_tools.h"

#include <algorithm>

idx
FastaReader::read_one_seq(Sequence &seq)
{
    seq.clear();
    bool need_defline = true;
    while (++m_reader) {
        const OneDataLine& odl = m_reader.get_line();
        if (odl.empty()) continue;
        const int c = odl.front();
        if (c == '>' || c == '@') {
            if (need_defline) {
                x_parse_defline(odl, seq.header());
                need_defline = false;
                continue;
            } else {
                m_reader.unget_line();
                break;
            }
        } else if (c == '+') {
            if (!++m_reader) mc_error << "quality score line is missing at around line " << m_reader.line_number() << eolog;
            break;
        } else if (is_comment_line(odl)) {
            continue;
        } else if (need_defline) {
			mc_error << "header line is missing at around line " << m_reader.line_number() << eolog;
        }

        x_parse_data_line(odl, seq.sequence());
    }

    if (seq.sequence().empty() && !seq.header().empty())
		mc_error << "near line " << m_reader.line_number() << ", sequence data is missing." << eolog;

    if (seq.header().empty() && seq.sequence().empty()) return -1;
    return seq.size();
}

void
FastaReader::x_parse_data_line(const OneDataLine& line, str_t& seq)
{
    x_check_data_line(line);
    const idx size = line.size();
    idx curr_pos = seq.size();
    idx pos = 0;
    for (pos = 0; pos < size; ++pos) {
        const int c = line[pos];
        if (c == ';') break;
        if (is_nucl(c) || c == '-') {
            seq.push_back(c);
        } else if (!isspace(c)) {
			mc_error << "There are invalid residue(s) around position "
					 << pos + 1
					 << " of line "
					 << m_reader.line_number()
					 << ": " 
					 << (char)c
					 << eolog;
        }
    }
}

void
FastaReader::x_parse_defline(const OneDataLine &line, str_t &header)
{
    header.clear();
    header.insert(header.end(), line.begin() + 1, line.end());
    if (header.empty()) {
		mc_error << "A sequence is given with an empty header around line "
				 << m_reader.line_number()
				 << eolog;
    }
}

void FastaReader::x_check_data_line(const OneDataLine &line)
{
    idx good = 0, bad = 0, size = line.size();
    idx ambig_nucl = 0;
    for (idx pos = 0; pos < size; ++pos) {
        const unsigned  char c = line[pos];
        if (is_alpha(c) || c == '*') {
            ++good;
            if (is_nucl(c) && is_ambig_nucl(c)) ++ambig_nucl;
        } else if (c == '-') {
            ++good;
        } else if (isspace(c) || (c >= '0' && c <= '9')) {

        } else if (c == ';') {
            break;
        } else {
            ++bad;
        }
    }

    if (bad >= good / 3 && (size > 3 || good == 0 || bad > good)){
		mc_error << "Near line "
				 << m_reader.line_number()
				 << ", there is a line that does not look like plausible data,"
				 << " but it's not marked as defline or comment line."
				 << eolog;
    }
}
