//
// Created by sysu on 2/26/17.
//

#include "buffer_line_reader.h"
#include <errno.h>
#include <cstring>

#include "aux_tools.h"

using namespace std;

BufferLineReader::BufferLineReader(const char *path)
{
    if (strcmp(path, "-") == 0) {
        m_file = gzdopen(fileno(stdin), "rb");
		if (!m_file) mc_error << "failed to open file '" << path << "': out of memory." << eolog;
    } else if ((m_file = gzopen(path, "rb")) == 0) {
		const char* err = errno ? strerror(errno) : "out of memory";
		mc_error << "failed to open file '" << path << "': " << err << eolog;
    }

    m_buf = new char[kBufferSize];
    m_done = false;
    m_unget_line = false;
    m_line_number = 0;
    x_read_buffer();
}

bool
BufferLineReader::eof() const
{
    return m_done && (m_cur >= m_buf_sz);
}

bool
BufferLineReader::operator++()
{
    ++m_line_number;
    if (m_unget_line) {
        m_unget_line = false;
        return true;
    }

#define buffer_read_ret (!(m_done && m_line.empty()))

    m_line.clear();
    const idx start = m_cur;
    const idx end = m_buf_sz;
    for (idx p = start; p < end; ++p) {
        const int c = m_buf[p];
        if (c == '\n') {
            m_line.insert(m_line.end(), m_buf + start, m_buf + p);
            m_cur = ++p;
            if (p == end) x_read_buffer();
            return buffer_read_ret;
        } else if (c == '\r') {
            m_line.insert(m_line.end(), m_buf + start, m_buf + p);
            m_cur = ++p;
            if (p == end) {
                if (x_read_buffer()) {
                    p = m_cur;
                    if (m_buf[p] == '\n') m_cur = p + 1;
                }
                return buffer_read_ret;
            }
            if (m_buf[p] != '\n') return buffer_read_ret;
            m_cur = ++p;
            if (p == end) x_read_buffer();
            return buffer_read_ret;
        }
    }

    x_load_long();
    return buffer_read_ret;
}

void
BufferLineReader::x_load_long()
{
    idx start = m_cur;
    idx end = m_buf_sz;
    m_line.insert(m_line.end(), m_buf + start, m_buf + end);
    while (x_read_buffer()) {
        start = m_cur;
        end = m_buf_sz;
        for (idx p = start; p < end; ++p) {
            const int c = m_buf[p];
            if (c == '\r' || c == '\n') {
                m_line.insert(m_line.end(), m_buf + start, m_buf + p);
                if (++p == end) {
                    if (x_read_buffer()) {
                        p = m_cur;
                        end = m_buf_sz;
                        if (p < end && c == '\r' && m_buf[p] == '\n') {
                            ++p;
                            m_cur = p;
                        }
                    }
                } else {
                    if (c == '\r' && m_buf[p] == '\n') {
                        if (++p == end) {
                            x_read_buffer();
                            p = m_cur;
                        }
                    }
                    m_cur = p;
                }
                return;
            }
        }
        m_line.insert(m_line.end(), m_buf + start, m_buf + end);
    }
}

bool
BufferLineReader::x_read_buffer()
{
    int r = gzread(m_file, m_buf, kBufferSize);
    if (r < 0) {
        int errnum;
        const char* errmsg = gzerror(m_file, &errnum);
		const char* why = (Z_ERRNO == errnum) ? strerror(errno) : errmsg;
		mc_error << "read error: " << why << eolog;
    }
    m_buf_sz = r;
    m_cur = 0;
    if (m_buf_sz == 0) {
        m_done = true;
        return false;
    } else if (m_buf_sz < kBufferSize) {
        m_done = true;
    }
    return true;
}

BufferLineReader::~BufferLineReader()
{
    gzclose(m_file);
    delete[] m_buf;
}

