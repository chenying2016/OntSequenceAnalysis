#include "align_tags.h"

bool
get_cns_tags(const std::string& qaln, 
			 const std::string& taln, 
			 int toff, 
			 double weight,
			 std::vector<AlignTag>& tags)
{
    AlignTag tag; tag.weight = weight;
    int jj = 0;
    int j = toff - 1;
    int p_j = -1;
    int p_jj = 0;
    char p_q_base = '-';
	
	for (std::size_t i = 0; i != qaln.size(); ++i) {
		if (qaln[i] != GAP_CHAR) {
			++jj;
		}
		if (taln[i] != GAP_CHAR) {
			jj = 0;
		}
		if (jj >= U8_MAX || p_jj >= U8_MAX) return 0;
		p_jj = jj;
	}

	jj = 0;
	p_jj = 0;
    for (std::size_t i = 0; i != qaln.size(); ++i) {
        if (qaln[i] != GAP_CHAR) ++jj;
        if (taln[i] != GAP_CHAR) {
			++j;
			jj = 0;
        }
        
        tag.t_pos = j;
        tag.p_t_pos = p_j; 
        tag.delta = (u8)jj;
        tag.p_delta = (u8)p_jj;
        tag.q_base = qaln[i];
        tag.p_q_base = p_q_base;

        p_j = j;
        p_jj = jj;
        p_q_base = qaln[i];

        tags.push_back(tag);
    }

    return 1;
}
