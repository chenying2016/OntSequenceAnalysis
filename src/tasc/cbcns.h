#ifndef CBCNS_H
#define CBCNS_H

#include <string>
#include <vector>

#include "../common/output_stream.h"
#include "align_tags.h"
#include "soa.h"
#include "cns_aux.h"

class CntBaseConsensus
{
public:
	CntBaseConsensus():
		li_alloc( sizeof(LinkInfo) ),
		dci_alloc( sizeof(DeltaCovInfo) ) {
		tags.reserve(1000000);
		backbone.reserve(MAX_SEQ_SIZE);
		coverage.reserve(MAX_SEQ_SIZE);
		}
		
	void clear() {
		tags.clear();
		backbone.clear();
		coverage.clear();
		li_alloc.Clear();
		dci_alloc.Clear();
	}
	
	void add_one_align(const std::string& qaln, 
					   const std::string& taln, 
					   int toff, 
					   double weight,
					   bool need_normalise_gaps);
	
	bool consensus(const int min_cov, 
				   const int min_size, 
				   const int template_id,
				   const int template_size, 
				   std::vector<std::pair<int, int> >& cns_intvs,
				   OutputStream& out);
	
	bool consensus(const int min_cov,
				   const int min_size,
				   const std::vector<char>& raw_read,
				   std::string& cns_seq);

private:
	std::vector<AlignTag>		tags;
	std::vector<BackboneItem>	backbone;
	std::vector<int>			coverage;
	LinkInfoAllocator			li_alloc;
	DeltaCovInfoAllocator		dci_alloc;
	std::string					nqaln;
	std::string					ntaln;
};

#endif // CBCNS_H
