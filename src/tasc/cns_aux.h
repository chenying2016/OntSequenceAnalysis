#ifndef CNS_AUX_H
#define CNS_AUX_H

#include "../common/defs.h"
#include "../common/mc_log.h"
#include "align_tags.h"
#include "soa.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

struct LinkInfo
{
	double weight;
    int p_t_pos;
    u8 p_delta;
    char p_q_base;
    int link_count;

    bool operator==(const LinkInfo& r) {
        return p_t_pos == r.p_t_pos
               &&
               p_delta == r.p_delta
               &&
               p_q_base == r.p_q_base;
    }
};

typedef FixedSizeObjectAllocator LinkInfoAllocator;

struct BaseLinks
{
    int n_link;
    int coverage;
    LinkInfo* plinks;
    int best_p_t_pos;
    u8 best_p_delta;
    u8 best_p_q_base;
    double score;
};

std::ostream& operator<<(std::ostream& out, const BaseLinks& lnk);

struct DeltaCovInfo
{
    BaseLinks links[5];
};

typedef FixedSizeObjectAllocator DeltaCovInfoAllocator;

struct BackboneItem
{
    int n_delta;
    DeltaCovInfo* delta;
};

std::ostream& operator<<(std::ostream& out, const BackboneItem& it);

struct BackboneItemCleaner
{
    void operator()(BackboneItem& bitem) {
        bitem.n_delta = 0;
        bitem.delta = 0;
    }
};

struct AlignTagPLinkEq
{
    bool operator()(const AlignTag& a, const AlignTag& b) {
        return a.p_t_pos == b.p_t_pos
               &&
               a.p_delta == b.p_delta
               &&
               a.p_q_base == b.p_q_base;
    }
};

void build_backbone(std::vector<AlignTag>& tags, 
					const int template_size, 
					DeltaCovInfoAllocator& dci_alloc,
					LinkInfoAllocator& li_alloc,
					std::vector<BackboneItem>& backbone,
					std::vector<int>& coverage);

void
consensus_backbone_segment(std::vector<BackboneItem>& backbone,
						   int from,
						   int to,
						   std::vector<int>& coverage,
						   std::string& cns_seq);

void
consensus_backbone_segment(std::vector<BackboneItem>& backbone,
						   int from,
						   int to,
						   std::vector<int>& coverage,
						   std::string& cns_seq,
						   int& cns_from,
						   int& cns_to);

#endif // CNS_AUX_H
