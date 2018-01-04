#include "mc_ostream.h"

namespace ns_mc_ostream
{

static const char digits[] = "9876543210123456789";
static const char* zero = digits + 9;

static const char digitsHex[] = "0123456789abcdef";

// efficient integer to string conversion, by Matthew Wilson
template <typename T>
int convert(char buf[], T value)
{
    T i = value;
    char* p = buf;

    do {
        int lsd = i % 10;
        i /= 10;
        *p++ = zero[lsd];
    } while (i);
    
    if (value < 0) {
        *p++ = '-';
    }
    *p = '\0';
    std::reverse(buf, p);    

    return static_cast<int>(p - buf);
}

template <typename T>
void
MCOstream::formatInteger(T v)
{
    char p[64];
    int len = convert(p, v);
    buffer.push_back(p, len);
}

MCOstream&
MCOstream::operator<<(short v)
{
	(*this) << static_cast<int>(v);
	return *this;
}

MCOstream&
MCOstream::operator<<(unsigned short v)
{
	(*this) << static_cast<unsigned int>(v);
	return *this;
}

MCOstream&
MCOstream::operator<<(int v)
{
	formatInteger(v);
	return *this;
}

MCOstream&
MCOstream::operator<<(unsigned int v)
{
	formatInteger(v);
	return *this;
}

MCOstream&
MCOstream::operator<<(long v)
{
    formatInteger(v);
    return *this;
}

MCOstream&
MCOstream::operator<<(unsigned long v)
{
    formatInteger(v);
    return *this;
}

MCOstream&
MCOstream::operator<<(long long v)
{
    formatInteger(v);
    return *this;
}

MCOstream&
MCOstream::operator<<(unsigned long long v)
{
    formatInteger(v);
    return *this;
}

MCOstream&
MCOstream::operator<<(double v)
{
	char p[64];
    int len = snprintf(p, 64, "%.12g", v);
    buffer.push_back(p, len);
    return *this;
}

};
