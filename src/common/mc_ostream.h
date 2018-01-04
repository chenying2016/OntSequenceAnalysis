#ifndef MC_OSTREAM_H
#define MC_OSTREAM_H

#include <string>

#include <cstring>

#include "pod_darr.h"

namespace ns_mc_ostream
{

class MCOstream
{
public:
	MCOstream(): buffer(1<<20) {}
	
	void append(const char* p, int len) {
		buffer.push_back(p, len);
	}
	
	void clear() {
		buffer.clear();
	}
	
	MCOstream& operator<<(bool v) {
		buffer.push_back(v ? '1' : '0');
		return *this;
	}
	
	MCOstream& operator<<(short);
	MCOstream& operator<<(unsigned short);
	MCOstream& operator<<(int);
	MCOstream& operator<<(unsigned int);
	MCOstream& operator<<(long);
	MCOstream& operator<<(unsigned long);
	MCOstream& operator<<(long long);
	MCOstream& operator<<(unsigned long long);
	
	MCOstream& operator<<(float v) {
		(*this) << static_cast<double>(v);
		return *this;
	}
	
	MCOstream& operator<<(double);
	
	MCOstream& operator<<(char v) {
		buffer.push_back(v);
		return *this;
	}
	
	MCOstream& operator<<(const char* v) {
		int n = strlen(v);
		buffer.push_back(v, n);
		return *this;
	}
	
	MCOstream& operator<<(const std::string& v) {
		buffer.push_back(v.c_str(), v.size());
		return *this;
	}
	
	char* data() {
		return buffer.data();
	}
	
	int size() {
		return buffer.size();
	}
	
private:
	template <typename T>
	void formatInteger(T v);
	
	PODArray<char> buffer;
};

}

#endif // MC_OSTREAM_H
