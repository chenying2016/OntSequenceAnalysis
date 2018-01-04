#include "overlaps_partition.h"

#include <cstdlib>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../common/aux_tools.h"
#include "../common/output_stream.h"
#include "../common/gapped_candidate.h"
#include "../common/pdb_aux.h"
#include "../common/record_reader.h"
#include "../common/smart_assert.h"

using namespace std;

void
make_partition_file_name(const char* prefix, const int part, string& name)
{
    name = prefix;
    ostringstream os;
    os << ".p" << part;
    name += os.str();
}

static int 
load_candidates(const char* can_input, vector<GappedCandidate>& canv)
{
    canv.clear();
    dsopen(ifstream, in, can_input, ios::in);
    GappedCandidate can;
    int max_id = -1;
    while (in >> can) {
        max_id = max(max_id, can.qid);
        max_id = max(max_id, can.sid);
        canv.push_back(can);
    }
    sclose(in);
    return max_id + 1;
}

void
init_parition_writers(const char* prefix, ofstream files[], OutputStream* pws[], const int spid, const int epid)
{
    int n = epid - spid;
    string name;
    for (int i = 0; i < n; ++i) {
        make_partition_file_name(prefix, spid + i, name);
        sopen(files[i], name.c_str(), ios::out);
        pws[i] = new OutputStream(files[i]);
    }
}

void
destroy_partition_writers(ofstream files[], OutputStream* pws[], const int n)
{
    for (int i = 0; i < n; ++i) {
        delete pws[i];
        sclose(files[i]);
    }
}

void
fix_candidate_offsets(const GappedCandidate& c, GappedCandidate& nc, const bool query_is_target)
{
	nc.score = c.score;
	
	if (query_is_target) {
		nc.qid = c.sid;
		nc.qsize = c.ssize;
		nc.sid = c.qid;
		nc.ssize = c.qsize;
		
		if (c.qdir == REV) {
			nc.qdir = REVERSE_STRAND(c.sdir);
			nc.qbeg = c.ssize - c.send;
			nc.qend = c.ssize - c.sbeg;
			nc.qoff = c.ssize - c.soff;
			nc.sdir = FWD;
			nc.sbeg = c.qsize - c.qend;
			nc.send = c.qsize - c.qbeg;
			nc.soff = c.qsize - c.qoff;
		} else {
			r_assert(c.qdir == FWD);
			nc.qdir = c.sdir;
			nc.qbeg = c.sbeg;
			nc.qend = c.send;
			nc.qoff = c.soff;
			nc.sdir = c.qdir;
			nc.sbeg = c.qbeg;
			nc.send = c.qend;
			nc.soff = c.qoff;
		}
	} else {
		nc.qid = c.qid;
		nc.qsize = c.qsize;
		nc.sid = c.sid;
		nc.ssize = c.ssize;
		
		if (c.sdir == REV) {
			nc.qdir = REVERSE_STRAND(c.qdir);
			nc.qbeg = c.qsize - c.qend;
			nc.qend = c.qsize - c.qbeg;
			nc.qoff = c.qsize - c.qoff;
			nc.sdir = FWD;
			nc.sbeg = c.ssize - c.send;
			nc.send = c.ssize - c.sbeg;
			nc.soff = c.ssize - c.soff;
		} else {
			r_assert(c.sdir == FWD);
			nc = c;
		}
	}
}

int
partition_candidates(const CnsOptions& options)
{
	const int batch_size = options.batch_size;
	const int min_read_size = options.min_size;
	const char* can_input = options.candidates;
    const int num_reads = load_num_reads(options.wrk_dir);
    const int num_batches = (num_reads + batch_size - 1) / batch_size;
    const int NumDumped = 10;
    ofstream output[NumDumped];
    OutputStream* pws[NumDumped];
	GappedCandidate can, ncan;

    for (int i = 0; i < num_batches; i += NumDumped) {
        const int sfid = i;
        const int efid = min(sfid + NumDumped, num_batches);
        const int nfid = efid - sfid;
        init_parition_writers(can_input, output, pws, sfid, efid);
        const int Lid = batch_size * sfid;
        const int Rid = batch_size * efid;
		
		RecordReader creader(can_input);
		while (creader.read_one_record(can)) {
            if (can.qsize < min_read_size || can.ssize < min_read_size) continue;
            if (can.qid >= Lid && can.qid < Rid) {
                int fid = (can.qid - Lid) / batch_size;
                r_assert(fid < nfid)(fid)(nfid)(can.qid)(Lid);
				fix_candidate_offsets(can, ncan, 1);
                (*pws[fid]) << ncan;
            } 
			if (can.sid >= Lid && can.sid < Rid) {
                int fid = (can.sid - Lid) / batch_size;
                r_assert(fid < nfid);
				fix_candidate_offsets(can, ncan, 0);
                (*pws[fid]) << ncan;
            }
		}

        destroy_partition_writers(output, pws, nfid);
    }

    return num_batches;
}
