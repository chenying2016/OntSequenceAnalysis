#include "consensus_aux.h"

#include <algorithm>

using namespace std;

void
build_cns_thrd_data(PackedDB* _reads,
					std::vector<GappedCandidate>* _candidates,
					std::size_t* _next_can_id,
					pthread_mutex_t* _can_id_lock,
					CnsOptions* _options,
					std::ostream* _out,
					pthread_mutex_t* _out_lock,
					std::ostream* _raw_out,
					pthread_mutex_t* _raw_out_lock,
					ConsensusData** ppctd)
{
	const int num_threads = _options->num_threads;
	for (int i = 0; i < num_threads; ++i) {
		ppctd[i] = new ConsensusData(_reads,
									 _candidates,
									 _next_can_id,
									 _can_id_lock,
									 _options,
									 _out,
									 _out_lock,
									 _raw_out,
									 _raw_out_lock);
	}
}

void
destroy_cns_thrd_data(ConsensusData** ppctd, const int num_threads)
{
	for (int i = 0; i < num_threads; ++i) {
		delete ppctd[i];
	}
}
