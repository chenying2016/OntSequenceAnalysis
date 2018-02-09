#ifndef CNS_SEQ_H
#define CNS_SEQ_H

#include "../klib/kstring.h"

typedef struct {
	int id;
	int left, right;
	int org_seq_size;
	kstring_t cns_seq;
} CnsSeq;

#define DUMP_CNS_SEQ(output_func, out, cns) \
	do { \
		const int __cns_read_size = (int)kstr_size((cns).cns_seq) - 1; \
		output_func(out, ">%d_%d_%d_%d_%d\n", (cns).id, (cns).org_seq_size, (cns).left, (cns).right, __cns_read_size); \
		DUMP_KSTRING(output_func, out, (cns).cns_seq); \
		output_func(out, "\n"); \
	} while (0)
   
#define clear_CnsSeq(cns) (kstr_clear((cns).cns_seq))
#define new_CnsSeq(cns) CnsSeq cns; kstr_init((cns).cns_seq)
#define free_CnsSeq(cns) free_kstring((cns).cns_seq)
#define init_CnsSeq(cns) kstr_init((cns).cns_seq)

#endif // CNS_SEQ_H
