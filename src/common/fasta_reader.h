//
// Created by sysu on 2/26/17.
//

#ifndef UNTITLED_FASTA_READER_H
#define UNTITLED_FASTA_READER_H

#include "buffer_line_reader.h"
#include "sequence.h"

class FastaReader
{
public:
    typedef Sequence::str_t         str_t;
    typedef BufferLineReader::OneDataLine   OneDataLine;

public:
    FastaReader(const char* path) : m_reader(path) { m_encode_table = get_dna_encode_table(); }
    idx read_one_seq(Sequence& seq);

private:
    void x_parse_defline(const OneDataLine& line, str_t& header);
    void x_parse_data_line(const OneDataLine& line, str_t& seq);
    void x_check_data_line(const OneDataLine& line);

    bool is_header_line(const OneDataLine& line) {
        return line.size() && (line[0] == '>' || line[0] == '@');
    }

    bool is_nucl(const unsigned char ch) {
        int r = m_encode_table[ch];
        return r < 16;
    }

    bool is_ambig_nucl(const unsigned char ch) {
        int r = m_encode_table[ch];
        return (r < 16) && (r > 3);
    }

    bool is_upper_case_letter(const char ch) {
        return ch >= 'A' && ch <= 'Z';
    }

    bool is_lower_case_letter(const char ch) {
        return ch >= 'a' && ch <= 'z';
    }

    bool is_alpha(const unsigned char ch) {
        return is_upper_case_letter(ch) || is_lower_case_letter(ch);
    }

    bool is_comment_line(const OneDataLine& line) {
        return line[0] == '#' || line[0] == '!';
    }

private:
    BufferLineReader    m_reader;
    const u8*           m_encode_table;
};

#endif //UNTITLED_FASTA_READER_H
