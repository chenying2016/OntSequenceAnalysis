//
// Created by sysu on 2/26/17.
//

#ifndef UNTITLED_BUFFER_LINE_READER_H
#define UNTITLED_BUFFER_LINE_READER_H

#include "defs.h"

#include <vector>
#include <zlib.h>

class BufferLineReader
{
public:
    typedef std::vector<char>   OneDataLine;

public:
    BufferLineReader(const char* path);
    ~BufferLineReader();

    OneDataLine& get_line() {
        return m_line;
    }

    void unget_line() {
        --m_line_number;
        m_unget_line = true;
    }

    idx line_number() const {
        return m_line_number;
    }

    bool eof() const;
    bool operator++();

private:
    void x_load_long();
    bool x_read_buffer();

private:
    static const idx    kBufferSize = 1024 * 1024 * 8;
    gzFile              m_file;
    char*               m_buf;
    idx                 m_cur;
    idx                 m_buf_sz;
    bool                m_done;
    OneDataLine         m_line;
    bool                m_unget_line;
    idx                 m_line_number;
};

#endif //UNTITLED_BUFFER_LINE_READER_H
