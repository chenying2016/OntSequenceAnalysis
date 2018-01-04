#ifndef ALIGN_TAGS_H
#define ALIGN_TAGS_H

#include "../common/defs.h"

#include <string>
#include <vector>

#define PB_WEIGHT 2
#define ONT_WEIGHT 1

struct AlignTag
{
	double weight;
	int t_pos;
    int p_t_pos;
    u8 delta;
    u8 p_delta;
    char q_base;
    char p_q_base;
};

struct AlignTagLessThan
{
    bool operator()(const AlignTag& a, const AlignTag& b) {
        if (a.t_pos != b.t_pos) return a.t_pos < b.t_pos;
        if (a.delta != b.delta) return a.delta < b.delta;
        if (a.q_base != b.q_base) return a.q_base < b.q_base;
        if (a.p_t_pos != b.p_t_pos) return a.p_t_pos < b.p_t_pos;
        if (a.p_delta != b.p_delta) return a.p_delta < b.p_delta;
        return a.p_q_base < b.p_q_base;
    }
};

bool
get_cns_tags(const std::string& qaln, 
			 const std::string& taln, 
			 int toff, 
			 double weight,
			 std::vector<AlignTag>& tags);

#endif // ALIGN_TAGS_H
