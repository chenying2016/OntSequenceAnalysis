#include "cns_aux.h"

using namespace std;

std::ostream& operator<<(std::ostream& out, const BaseLinks& lnk)
{
	out << "n_link = " << lnk.n_link
		<< ", coverage = " << lnk.coverage
		<< ", best_p_t_pos = " << lnk.best_p_t_pos 
		<< ", best_p_delta = " << (int)lnk.best_p_delta
		<< ", best_p_q_base = " << (int)lnk.best_p_q_base
		<< ", score = " << lnk.score
		<< '\n';
	return out;
}

std::ostream& operator<<(std::ostream& out, const BackboneItem& it)
{
	out << "n_delta = " << it.n_delta << '\n'
		<< "\tlinks[0]: " << it.delta->links[0]
		<< "\tlinks[1]: " << it.delta->links[1]
		<< "\tlinks[2]: " << it.delta->links[2]
		<< "\tlinks[3]: " << it.delta->links[3]
		<< "\tlinks[4]: " << it.delta->links[4]
		<< '\n';
	return out;
}

inline u8
encode_dna_base(const char c)
{
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case '-': return 4;
    }
    mc_error << "illegal dna base: '" << c << "'" << eolog;
    return 4; // we will never come here, but the compiler will complain if we omit it.
}

void
build_base_links(AlignTag* tags, const int ntag, BaseLinks& link, LinkInfoAllocator& li_alloc)
{
    int n_link = 0;
    int i = 0;
    AlignTagPLinkEq tag_plink_eq;

    while (i < ntag) {
        int j = i + 1;
        while (j < ntag && tag_plink_eq(tags[i], tags[j])) ++j;
        ++n_link;
        i = j;
    }
    
    link.plinks = (LinkInfo*)li_alloc.Allocate(n_link);
    link.n_link = n_link;
    link.coverage = ntag;

    n_link = 0;
    i = 0;
    while (i < ntag) {
        int j = i + 1;
        while (j < ntag && tag_plink_eq(tags[i], tags[j])) ++j;
        LinkInfo& linfo = link.plinks[n_link];
        linfo.p_t_pos = tags[i].p_t_pos;
        linfo.p_delta = tags[i].p_delta;
        linfo.p_q_base = tags[i].p_q_base;
        linfo.link_count = j - i;
		linfo.weight = 0;
		for (int k = i; k < j; ++k) {
			linfo.weight += tags[k].weight;
		}

        ++n_link;
        i = j;
    }

    r_assert(n_link == link.n_link);
}

void
build_delta_links(AlignTag* tags, const int ntag, DeltaCovInfo& dci, LinkInfoAllocator& li_alloc)
{
    for (int i = 0; i < 5; ++i) {
        dci.links[i].n_link = 0;
        dci.links[i].coverage = 0;
        dci.links[i].best_p_t_pos = -1;
        dci.links[i].best_p_delta = U8_MAX;
        dci.links[i].best_p_q_base = '.';
        dci.links[i].score = 0;
    }

    int i = 0;
    while (i < ntag) {
        int j = i;
        while (j < ntag && tags[i].q_base == tags[j].q_base) ++j;
        u8 c = encode_dna_base(tags[i].q_base);
        build_base_links(tags + i, j - i, dci.links[c], li_alloc);
        i = j;
    }
}

void
build_backbone_item(AlignTag* tags, 
					const int ntag, 
					BackboneItem& bitem, 
					DeltaCovInfoAllocator& dci_alloc, 
					LinkInfoAllocator& li_alloc,
				    vector<int>& coverage)
{
    bitem.n_delta = tags[ntag - 1].delta + 1;
    bitem.delta = (DeltaCovInfo*)dci_alloc.Allocate(bitem.n_delta);
    int i = 0;
    while (i < ntag) {
        int j = i + 1;
        while (j < ntag && tags[i].delta == tags[j].delta) ++j;
        build_delta_links(tags + i, j - i, bitem.delta[tags[i].delta], li_alloc);
		if (tags[i].delta == 0) {
			r_assert(tags[i].t_pos < coverage.size());
			coverage[ tags[i].t_pos ] = j - i;
		}
        i = j;
    }
}

void build_backbone(std::vector<AlignTag>& tags, 
					const int template_size, 
					DeltaCovInfoAllocator& dci_alloc,
					LinkInfoAllocator& li_alloc,
					std::vector<BackboneItem>& backbone,
					std::vector<int>& coverage)
{
    backbone.resize(template_size);
    std::for_each(backbone.begin(), backbone.end(), BackboneItemCleaner());
    coverage.resize(template_size);
    fill(coverage.begin(), coverage.end(), 0);

    std::sort(tags.begin(), tags.end(), AlignTagLessThan());
    int i = 0, n = tags.size();
    while (i < n) {
        int j = i + 1;
        while (j < n && tags[i].t_pos == tags[j].t_pos) ++j;
		r_assert(tags[i].t_pos < template_size)(i)(tags[i].t_pos)(template_size);
			 
		r_assert(i < tags.size())(i)(tags.size());
		r_assert(j <= tags.size())(j)(tags.size());
        build_backbone_item(tags.data() + i, j - i, backbone[tags[i].t_pos], dci_alloc, li_alloc, coverage);
        i = j;
    }
}

void
consensus_backbone_segment(std::vector<BackboneItem>& backbone,
						   int from,
						   int to,
						   vector<int>& coverage,
						   std::string& cns_seq)
{
    int g_best_ck = 0;
    BaseLinks* g_best_aln_col = 0;
    int g_best_t_pos = 0;
    double g_best_score = -1.0;
	int g_best_q_base = -1;
    int kk;
    int ck;
    int best_i;
    int best_j;
    int best_b;
    int best_ck = -1;
    double score;
    double best_score;
    BaseLinks* aln_col;

    for (int i = from; i < to; ++i) {
        for (int j = 0; j < backbone[i].n_delta; ++j) {
            for (kk = 0; kk < 5; ++kk) {
                aln_col = backbone[i].delta[j].links + kk;
                if (aln_col->coverage) {
                    best_score = -1;
                    best_i = -1;
                    best_j = -1;
                    best_b = -1;
                    for (ck = 0; ck < aln_col->n_link; ++ck) {
                        int pi = aln_col->plinks[ck].p_t_pos;
                        int pj = aln_col->plinks[ck].p_delta;
                        int pkk = encode_dna_base(aln_col->plinks[ck].p_q_base);
						//if (j) score = 1.0 * aln_col->plinks[ck].link_count - 0.4 * 0.5 * coverage[i];
						//else score = 1.0 * aln_col->plinks[ck].link_count - 0.5 * coverage[i];
						if (j) score = 1.0 * aln_col->plinks[ck].weight - 0.4 * 0.5 * coverage[i];
						else score = 1.0 * aln_col->plinks[ck].weight - 0.4 * 0.5 * coverage[i];
                        if (aln_col->plinks[ck].p_t_pos != -1) {
                            score += backbone[pi].delta[pj].links[pkk].score;
                        }

                        if (score > best_score) {
                            best_score = score;
                            aln_col->best_p_t_pos = best_i = pi;
                            aln_col->best_p_delta = best_j = pj;
                            aln_col->best_p_q_base = best_b = pkk;
                            best_ck = ck;
                        }
                    }

                    aln_col->score = best_score;
                    if (best_score > g_best_score) {
                        g_best_score = best_score;
                        g_best_aln_col = aln_col;
                        g_best_ck = best_ck;
                        g_best_t_pos = i;
						g_best_q_base = kk;
                    }
                }
            }
        }
    }
    r_assert(g_best_score != -1);

    std::string t_cns_seq;
    char bb = '$';
    ck = g_best_q_base;
    int i = g_best_t_pos;
    int j;

    while (1) {
        bb = DNADecoder::decode(ck); 
        i = g_best_aln_col->best_p_t_pos;
        if (i == -1) break;
        j = g_best_aln_col->best_p_delta;
        ck = g_best_aln_col->best_p_q_base;
        g_best_aln_col = backbone[i].delta[j].links + ck;
        
        if (bb != '-') t_cns_seq += bb;
    }

    cns_seq.assign(t_cns_seq.rbegin(), t_cns_seq.rend());
}

void
consensus_backbone_segment(std::vector<BackboneItem>& backbone,
						   int from,
						   int to,
						   vector<int>& coverage,
						   std::string& cns_seq,
						   int& cns_from,
						   int& cns_to)
{
    int g_best_ck = 0;
    BaseLinks* g_best_aln_col = 0;
    int g_best_t_pos = 0;
    double g_best_score = -1.0;
	int g_best_q_base = -1;
    int kk;
    int ck;
    int best_i;
    int best_j;
    int best_b;
    int best_ck = -1;
    double score;
    double best_score;
    BaseLinks* aln_col;

    for (int i = from; i < to; ++i) {
        for (int j = 0; j < backbone[i].n_delta; ++j) {
            for (kk = 0; kk < 5; ++kk) {
                aln_col = backbone[i].delta[j].links + kk;
                if (aln_col->coverage) {
                    best_score = -1;
                    best_i = -1;
                    best_j = -1;
                    best_b = -1;

                    for (ck = 0; ck < aln_col->n_link; ++ck) {
                        int pi = aln_col->plinks[ck].p_t_pos;
                        int pj = aln_col->plinks[ck].p_delta;
                        int pkk = encode_dna_base(aln_col->plinks[ck].p_q_base);
						if (j) score = 1.0 * aln_col->plinks[ck].link_count - 0.4 * 0.5 * coverage[i];
						else score = 1.0 * aln_col->plinks[ck].link_count - 0.5 * coverage[i];
                        if (aln_col->plinks[ck].p_t_pos != -1) {
                            score += backbone[pi].delta[pj].links[pkk].score;
                        }

                        if (score > best_score) {
                            best_score = score;
                            aln_col->best_p_t_pos = best_i = pi;
                            aln_col->best_p_delta = best_j = pj;
                            aln_col->best_p_q_base = best_b = pkk;
                            best_ck = ck;
                        }
                    }

                    aln_col->score = best_score;
                    if (best_score > g_best_score) {
                        g_best_score = best_score;
                        g_best_aln_col = aln_col;
                        g_best_ck = best_ck;
                        g_best_t_pos = i;
						g_best_q_base = kk;
                    }
                }
            }
        }
    }
    r_assert(g_best_score != -1);

    std::string t_cns_seq;
    char bb = '$';
    ck = g_best_q_base;
    int i = g_best_t_pos;
	cns_to = i + 1;
    int j;

    while (1) {
		bb = ck;
        i = g_best_aln_col->best_p_t_pos;
        if (i == -1) break;
        j = g_best_aln_col->best_p_delta;
        ck = g_best_aln_col->best_p_q_base;
        g_best_aln_col = backbone[i].delta[j].links + ck;
		cns_from = i;
		
		if (bb != 4) {
			r_assert(bb >= 0 && bb < 4);
			
			t_cns_seq += bb;
		}
    }

    cns_seq.assign(t_cns_seq.rbegin(), t_cns_seq.rend());
}
