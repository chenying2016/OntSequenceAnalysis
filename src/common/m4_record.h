#ifndef M4_RECORD_H
#define M4_RECORD_H

#include "../common/ontcns_defs.h"
#include "../klib/kvec.h"

typedef struct {
	int qid;
    int qdir;
    idx qoff;
    idx qend;
    idx qext;
    idx qsize;
    int sid;
    int sdir;
    idx soff;
    idx send;
    idx sext;
    idx ssize;
    double ident_perc;
    int vscore;
} M4Record;

typedef kvec_t(M4Record) vec_m4;

#define DUMP_M4_RECORD(output_func, out, m) \
	output_func(out, "%d\t%d\t%.2f\t%d\t%d\t%lu\t%lu\t%lu\t%lu\t%d\t%lu\t%lu\t%lu\t%lu\n", \
				(m).qid, \
				(m).sid, \
				(m).ident_perc, \
				(m).vscore, \
				(m).qdir, \
				(m).qoff, \
				(m).qend, \
				(m).qext, \
				(m).qsize, \
				(m).sdir, \
				(m).soff, \
				(m).send, \
				(m).sext, \
				(m).ssize)
   
#define LOAD_M4_RECORD(input_func, in, m) \
	input_func(in, "%d%d%lf%d%d%lu%lu%lu%lu%d%lu%lu%lu%lu", \
			   &(m).qid, \
			   &(m).sid, \
			   &(m).ident_perc, \
			   &(m).vscore, \
			   &(m).qdir, \
			   &(m).qoff, \
			   &(m).qend, \
			   &(m).qext, \
			   &(m).qsize, \
			   &(m).sdir, \
			   &(m).soff, \
			   &(m).send, \
			   &(m).sext, \
			   &(m).ssize)

typedef enum {
	eM4_Contained,
	eM4_Contains,
	eM4_Overlap,
	eM4_None
} EM4OverlapType;

EM4OverlapType
detect_m4_type(const M4Record* m, const idx unmapped_boundary);

#endif // M4_RECORD_H
