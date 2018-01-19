#ifndef KSTRING_H
#define KSTRING_H

#include <stdlib.h>
#include <string.h>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#define KSTR_SIZE(str) ((str).l)
#define KSTR_DATA(str) ((str).s)
#define KSTR_INIT(str) do { (str).l = 0; (str).m = 0; (str).s = 0; } while(0)
#define KSTR_DESTROY(str) do { free((str).s); (str).l = 0; (str).m = 0; (str).s = 0; } while(0)
#define KSTR_CLEAR(str) do { (str).l = 0; } while(0)

kstring_t*
new_kstring();

kstring_t*
free_kstring(kstring_t* s);

void kstring_set(kstring_t* s, size_t pos, char c);

void kstring_reserve(kstring_t *s, size_t size);

void kstring_resize(kstring_t* s, size_t size);

int kputsn(const char *p, int l, kstring_t *s);

int kputs(const char *p, kstring_t *s);

int kputc(int c, kstring_t *s);

int kputw(int c, kstring_t *s);

int kputuw(unsigned c, kstring_t *s);

int kputl(long c, kstring_t *s);

int ksprintf(kstring_t *s, const char *fmt, ...);

#endif
