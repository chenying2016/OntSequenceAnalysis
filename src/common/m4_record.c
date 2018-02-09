#include "m4_record.h"

EM4OverlapType
detect_m4_type(const M4Record* m, const idx unmapped_boundary)
{
	if (m->qoff <= unmapped_boundary && m->qsize - m->qend <= unmapped_boundary) return eM4_Contained;
	if (m->soff <= unmapped_boundary && m->ssize - m->send <= unmapped_boundary) return eM4_Contains;
	if (m->qsize - m->qend <= unmapped_boundary && m->soff <= unmapped_boundary) return eM4_Overlap;
	if (m->ssize - m->send <= unmapped_boundary && m->qoff <= unmapped_boundary) return eM4_Overlap;
	return eM4_None;
}
