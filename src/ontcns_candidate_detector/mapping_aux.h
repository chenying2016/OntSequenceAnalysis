#ifndef MAPPING_AUX_H
#define MAPPING_AUX_H

#include "../common/aux_tools.h"
#include "../common/lookup_table.h"
#include "../common/packed_db.h"
#include "../common/output_stream.h"
#include "../word_finder/word_finder_aux.h"
#include "options.h"


struct MappingData
{
   int              thread_id;
   MapOptions*		options;
   PackedDB*        reads;
   int              reads_start_id;
   PackedDB*        reference;
   int              reference_start_id;
   LookupTable*     lktbl;
   OutputStream	    out;
   const int        chunk_size;
   int*             chunk_id;
   int              num_chunks;
   int              num_reads;
   pthread_mutex_t* chunk_lock;

   MappingData(int tid, MapOptions* mo, PackedDB* _reads,
           PackedDB* _reference, LookupTable* lt, 
           std::ostream& _out, pthread_mutex_t* out_lock,
		   int* _chunk_id, pthread_mutex_t* _chunk_lock)
       : thread_id(tid),
         options(mo),
         reads(_reads),
         reference(_reference),
         lktbl(lt),
         out(_out, out_lock),
		 chunk_size(500),
		 chunk_id(_chunk_id),
		 chunk_lock(_chunk_lock) {
         }
		 
	void reset_reads_info(int reads_sid, int reference_sid) {
		reads_start_id = reads_sid;
		reference_start_id = reference_sid;
		num_reads = (int)reads->num_seqs();
		num_chunks = (num_reads + chunk_size - 1) / chunk_size;
	}

   bool get_next_chunk(int& sid, int& eid) {
       bool r = true;
       pthread_mutex_lock(chunk_lock);
	   int cid = *chunk_id;
       if (cid >= num_chunks) {
           r = false;
       } else {
           sid = cid * chunk_size;
           eid = sid + chunk_size;
           eid = std::min(eid, num_reads);
           ++*chunk_id;
       }
       pthread_mutex_unlock(chunk_lock);
       return r;
   }
};

#endif // MAPPING_AUX_H
