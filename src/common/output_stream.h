#ifndef OUTPUT_STREAM_H
#define OUTPUT_STREAM_H

#include <cassert>
#include <cstring>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "aux_tools.h"
#include "mc_ostream.h"

namespace ns_mc_output_stream
{

using std::ofstream;
using std::string;
using std::ostream;
using ns_mc_ostream::MCOstream;
typedef pthread_mutex_t Mutex;

class OutputStream: public mc_noncopyable
{
    typedef OutputStream self;
public:
    OutputStream(ostream& file, Mutex* file_lock = NULL):
        ext_file_(&file),
        ext_file_lock_(file_lock),
        use_int_file_(0) {
		snew(buffer_, char, kBufferSize);
		cur_ = 0;
	}

    OutputStream(const char* path) {
        snew(buffer_, char, kBufferSize);
        cur_ = 0;
        sopen(int_file_, path, std::ios::out);
        ext_file_ = NULL;
        ext_file_lock_ = NULL;
        use_int_file_ = 1;
    }
	
	~OutputStream() {
		dump();
		if (use_int_file_) {
			sclose(int_file_);
		}
		sfree(buffer_);
	}
	
	self& operator<<(bool v) {
		os << v;
		receive_os_data();
		return *this;
	}
	
	self& operator<<(char v) {
		os << v;
		receive_os_data();
		return *this;
	}
	
	self& operator<<(short v) {
		os << v;
		receive_os_data();
		return *this;
	}
	
	self& operator<<(unsigned short v) {
		os << v;
		receive_os_data();
		return *this;
	}
	
	self& operator<<(int v) {
		os << v;
		receive_os_data();
		return *this;
	}
	
	self& operator<<(unsigned int v) {
		os << v;
		receive_os_data();
		return *this;
	}
	
	self& operator<<(long v) {
		os << v;
		receive_os_data();
		return *this;
	}
	
	self& operator<<(unsigned long v) {
		os << v;
		receive_os_data();
		return *this;
	}
	
	self& operator<<(long long v) {
		os << v;
		receive_os_data();
		return *this;
	}
	
	self& operator<<(unsigned long long v) {
		os << v;
		receive_os_data();
		return *this;
	}
	
	self& operator<<(float v) {
		os << v;
		receive_os_data();
		return *this;
	}
	
	self& operator<<(double v) {
		os << v;
		receive_os_data();
		return *this;
	}
	
	self& operator<<(const char* v) {
		os << v;
		receive_os_data();
		return *this;
	}
	
	self& operator<<(const string& v) {
		os << v;
		receive_os_data();
		return *this;
	}
	
	MCOstream& mcostream() {
		return os;
	}
	
	void absorp_os_data() {
		receive_os_data();
	}

private:
    void dump() {
        if (!cur_) return;
        if (use_int_file_) {
            sbwrite(int_file_, buffer_, cur_);
            cur_ = 0;
        } else {
            if (ext_file_lock_) pthread_mutex_lock(ext_file_lock_);
            sbwrite(*ext_file_, buffer_, cur_);
            cur_ = 0;
            if (ext_file_lock_) pthread_mutex_unlock(ext_file_lock_);
        }
    }

    int avail() {
        return kBufferSize - cur_;
    }

    void receive_os_data() {
        if (os.size() > avail()) dump();
		memcpy(buffer_ + cur_, os.data(), os.size());
		cur_ += os.size();
		os.clear();
    }
    
private:
	char* buffer_;
	int cur_;
	static const int kBufferSize = 8 * 1024 * 1024;
	MCOstream os;

    ofstream int_file_;
    ostream* ext_file_;
    Mutex* ext_file_lock_;
    bool use_int_file_;
};

} // ns_mc_output_stream

typedef ns_mc_ostream::MCOstream MCOstream;
typedef ns_mc_output_stream::OutputStream OutputStream;

#endif // OUTPUT_STREAM_H
