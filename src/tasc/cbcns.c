#include "cbcns.h"

#include "../common/oc_assert.h"

RESULT_WRITER_IMPL(CnsSeq, CNS_SEQ, static) 

CbCnsData*
new_CbCnsData()
{
	CbCnsData* data = (CbCnsData*)malloc(sizeof(CbCnsData));
	kv_init(data->tags);
	kv_init(data->backbone);
	kv_init(data->coverage);
	data->li_alloc = new_OcObjectAllocator(sizeof(LinkInfo));
	data->dci_alloc = new_OcObjectAllocator(sizeof(DeltaCovInfo));
	return data;
}

CbCnsData*
free_CbCnsData(CbCnsData* data)
{
	kv_destroy(data->tags);
	kv_destroy(data->backbone);
	kv_destroy(data->coverage);
	data->li_alloc = free_OcObjectAllocator(data->li_alloc);
	data->dci_alloc = free_OcObjectAllocator(data->dci_alloc);
	free(data);
	return 0;
}

void
clear_CbCnsData(CbCnsData* data)
{
	kv_clear(data->tags);
	kv_clear(data->backbone);
	kv_clear(data->coverage);
	clear_OcObjectAllocator(data->li_alloc);
	clear_OcObjectAllocator(data->dci_alloc);
}

void
add_one_align(CbCnsData* cns_data,
			  const char* qaln,
			  const char* taln,
			  const size_t aln_size,
			  kstring_t* target,
			  int toff,
			  int tend,
			  double weight)
{
	get_cns_tags(qaln,
				 taln,
				 aln_size,
				target,
				 toff,
				 tend,
				 weight,
				 &cns_data->tags);
}

void
consensus_broken(CbCnsData* cns_data,
				  const int min_cov,
				  const int min_size,
				  const int template_id,
				  const int template_size,
				  vec_intpair* cns_intvs,
				  CnsSeq* cns_seq,
				  ResultsWriter* out)
{
	build_backbone(kv_data(cns_data->tags), 
				   kv_size(cns_data->tags), 
				   template_size, 
				   cns_data->dci_alloc, 
				   cns_data->li_alloc,
				   &cns_data->backbone,
				   &cns_data->coverage);
	
	int i = 0;
	int* coverage = kv_data(cns_data->coverage);
	while (i < template_size) {
		while (i < template_size && coverage[i] < min_cov) ++i;
		int j = i + 1;
		while (j < template_size && coverage[j] >= min_cov) ++j;
		if (j - i >= min_size * 0.85) {
			clear_CnsSeq(*cns_seq);
			consensus_backbone_segment(kv_data(cns_data->backbone),
									   i,
									   j,
									   coverage,
									   &cns_seq->cns_seq,
									   NULL,
									   NULL);
			if (kstr_size(cns_seq->cns_seq) >= min_size) {
				cns_seq->id = template_id;
				cns_seq->left = i;
				cns_seq->right = j;
				cns_seq->org_seq_size = template_size;
				for (size_t k = 0; k != kstr_size(cns_seq->cns_seq); ++k) {
					int c = kstr_A(cns_seq->cns_seq, k);
					kstr_A(cns_seq->cns_seq, k) = DecodeDNA(c);
				}
				RW_DUMP_ONE_DATA(CNS_SEQ, out, cns_seq);
				IntPair ip;
				ip.first = i;
				ip.second = j;
				kv_push(IntPair, *cns_intvs, ip);
			}
		}
		i = j;
	}
}

typedef struct {
	int raw_from, raw_to;
	int cns_from, cns_to;
} CnsIntv;
typedef kvec_t(CnsIntv) vec_cns_intv;

void
consensus_unbroken(CbCnsData* cns_data,
			   const int min_cov,
			   const int min_size,
			   const char* raw_read,
			   const int template_size,
			   kstring_t* cns_seq)
{
	build_backbone(kv_data(cns_data->tags), 
				   kv_size(cns_data->tags), 
				   template_size, 
				   cns_data->dci_alloc, 
				   cns_data->li_alloc,
				   &cns_data->backbone,
				   &cns_data->coverage);
	int i = 0;
	int raw_from, raw_to, cns_from, cns_to;
	new_kstring(cns_frag);
	new_kstring(cns_all);
	new_kvec(vec_cns_intv, cns_intvs);
	CnsIntv intv;
	int* coverage = kv_data(cns_data->coverage);
	kstr_clear(*cns_seq);
	while (i < template_size) {
		while (i < template_size && coverage[i] < min_cov) ++i;
		int j = i + 1;
		while (j < template_size && coverage[j] >= min_cov) ++j;
		
		if (j - i >= min_size * 0.85) {
			consensus_backbone_segment(kv_data(cns_data->backbone),
									   i,
									   j,
									   coverage,
									   &cns_frag,
									   &raw_from,
									   &raw_to);
			if (kstr_size(cns_frag) >= min_size) {
				cns_from = kstr_size(cns_all);
				cns_to = cns_from + kstr_size(cns_frag);
				kputsn(kstr_str(cns_frag), kstr_size(cns_frag), &cns_all);
				intv.raw_from = raw_from;
				intv.raw_to = raw_to;
				intv.cns_from = cns_from;
				intv.cns_to = cns_to;
				kv_push(CnsIntv, cns_intvs, intv);
			}
		}
		i = j;
	}
	
	if (kv_size(cns_intvs) == 0) return;
	
	int last_raw_to = 0;
	int last_cns_to = 0;
	for (size_t k = 0; k != kv_size(cns_intvs); ++k) {
		raw_from = kv_A(cns_intvs, k).raw_from;
		raw_to = kv_A(cns_intvs, k).raw_to;
		cns_from = kv_A(cns_intvs, k).cns_from;
		cns_to = kv_A(cns_intvs, k).cns_to;
		oc_assert(cns_from == last_cns_to);
		if (raw_from > last_raw_to) {
			oc_assert(raw_from < template_size);
			kputsn(raw_read + last_raw_to, raw_from - last_raw_to, cns_seq);
		}
		kputsn(kstr_str(cns_all) + cns_from, cns_to - cns_from, cns_seq);
		last_raw_to = raw_to;
		last_cns_to = cns_to;
		oc_assert(last_raw_to <= template_size);
	}
	
	if (last_raw_to != template_size) {
		oc_assert(last_raw_to < template_size);
		kputsn(raw_read + last_raw_to, template_size - last_raw_to, cns_seq);
	}
	
	for (size_t k = 0; k != kstr_size(*cns_seq); ++k) {
		int c = kstr_A(*cns_seq, k);
		oc_assert(c >= 0 && c < 4);
		kstr_A(*cns_seq, k) = DecodeDNA(c);
	}
}
