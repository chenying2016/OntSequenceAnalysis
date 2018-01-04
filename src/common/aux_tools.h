#ifndef AUX_TOOLS_H
#define AUX_TOOLS_H

#include "mc_log.h"

/// file manipulator

#define sopen(file, path, mode) \
do { \
    (file).open(path, mode); \
	if (!(file)) mc_error << "failed to open file " << path << " with mode " << #mode << eolog; \
} while (0)

#define dsopen(type, file, path, mode) type file; sopen(file, path, mode);

#define sclose(file) (file).close();

#define sbwrite(file, buf, size) \
do { \
    std::streamsize __sw__size = (size); \
    std::streambuf* __sw__sb = (file).rdbuf(); \
    const char* __sw__buf = (const char*)(buf); \
    std::streamsize __sw__ws = __sw__sb->sputn(__sw__buf, size); \
    if (__sw__ws != __sw__size) mc_error << "file write error!" << eolog; \
} while (0)

#define sbread(file, buf, size) \
do { \
    std::streamsize __sr__size = (size); \
    std::streambuf* __sr__sb = (file).rdbuf(); \
    char* __sr__buf = (char*)(buf); \
    std::streamsize __sr__rs = __sr__sb->sgetn(__sr__buf, size); \
    if (__sr__rs != __sr__size) mc_error << "file read error!" << eolog; \
} while (0)

/// memory alloc/dealloc
#include <new>
#define snew(arr, type, count) \
do { \
    size_t snew_cnt = static_cast<size_t>(count); \
    try { \
        (arr) = new type[snew_cnt]; \
    } catch (std::bad_alloc& snew_e) { \
		mc_error << snew_e.what() << eolog; \
    } \
} while(0)

#define sznew(arr, type, count) \
do { \
    size_t zsnew_cnt = static_cast<size_t>(count); \
    snew(arr, type, count); \
    size_t zsnew_s = zsnew_cnt * sizeof(type); \
    memset(static_cast<void*>(arr), 0, zsnew_s); \
} while(0)

#define dsnew(arr, type, count) type * arr = NULL; snew(arr, type, count)

#define dsznew(arr, type, count) type * arr = NULL; sznew(arr, type, count)

#define sfree(arr) delete[] arr

/// aux functions
inline bool
is_power_of_2(i64 n)
{
	if (n <= 0) return false;
	return (n & (n - 1)) == 0;
}

#include <sys/stat.h>
inline idx
file_size(const char* path)
{
    struct stat info;
    int r = stat(path, &info);
    if (r) {
        perror("stat");
    }
    return static_cast<idx>(info.st_size);
}

#include <algorithm>
// efficient integer to string conversion, by Matthew Wilson
template <typename T>
int int2string(char buf[], T value)
{
	const char* digits = "9876543210123456789";
	const char* zero = digits + 9;
	T i = value;
	char* p = buf;
	
	do {
		int lsd = i % 10;
		i /= 10;
		*p++ = zero[lsd];
	} while (i);
	
	if (value < 0) *p++ = '-';
	*p = '\0';
	std::reverse(buf, p);
	return static_cast<int>(p - buf);
}

struct mc_noncopyable
{
public:
    mc_noncopyable() {}

private:
    mc_noncopyable(const mc_noncopyable&);
    mc_noncopyable& operator=(const mc_noncopyable&);
};

inline void print_size(idx size)
{
	const idx one = 1;
	const idx GB = one << 30;
	const idx MB = one << 20;
	const idx KB = one << 10;
	
	if (size >= GB) {
		std::cout << (size >> 30) << "GB";
	} else if (size >= MB) {
		std::cout << (size >> 20) << "MB";
	} else if (size >= KB) {
		std::cout << (size >> 10) << "KB";
	} else {
		std::cout << size << "B";
	}
}

#endif // AUX_TOOLS_H
