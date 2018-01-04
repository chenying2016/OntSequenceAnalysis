#include "../common/aux_tools.h"
#include "../common/packed_db.h"
#include "../common/fasta_reader.h"
#include "../common/mc_log.h"
#include "../common/pdb_aux.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

static const idx kVolSize = DEFAULT_VOLUME_SIZE;

int
pack_one_file(const char* path, 
			  PackedDB& reads,
			  const char* wrk_dir,
			  const int platform,
			  int& read_id,
			  int& vid)
{
	int num_reads = 0;
	idx num_bps = 0;
	idx cvs = reads.size();
	FastaReader reader(path);
	Sequence read;
	while (1) {
		idx s = reader.read_one_seq(read);
		if (s == -1) break;
		++num_reads;
		num_bps += s;
		cvs += s;
		check_and_rename_ontcns_header(read, read_id++);
		reads.add_one_raw_seq(read, platform);
		if (cvs >= kVolSize) {
			string vname;
			make_volume_name(wrk_dir, vid++, vname);
			reads.dump_pac(vname.c_str());
			reads.clear();
            cvs = 0;
		}
	}
	
	mc_log << "packed " << path << ": " 
		   << num_reads << " reads (" << num_bps << " bps)"
		   << eolog;
	return num_reads;
}

int
pack_one_list(const char* list_path,
			  const int platform,
			  PackedDB& reads,
			  int& read_id,
			  int& vid,
			  const char* wrk_dir)
{
	dsopen(ifstream, in, list_path, ios::in);
	string line;
	int num_reads = 0;
	while (getline(in, line)) {
		num_reads += pack_one_file(line.c_str(), reads, wrk_dir, platform, read_id, vid); 
	}
	sclose(in);
	return num_reads;
}

void 
print_usage(const char* prog)
{
	cerr << "USAGE:\n"
		 << prog << ' ' << "wrk-dir file_list platform [file_list platform]" << endl;
}

int main(int argc, char* argv[])
{
	if (argc < 4) {
		print_usage(argv[0]);
		return 1;
	}
	
	PackedDB reads;
	reads.alloc(kVolSize + MAX_SEQ_SIZE);
	const char* wrk_dir = argv[1];
	int c = 2;
	int vid = 0;
	int num_reads = 0;
	int read_id = 0;
	while (c < argc) {
		const char* file_list = argv[c];
		int platform = atoi(argv[c + 1]);
		if (platform != TECH_PACBIO && platform != TECH_NANOPORE) {
			cerr << "invalid platform option: " << platform << endl;
			return 1;
		}
		num_reads += pack_one_list(file_list, platform, reads, read_id, vid, wrk_dir);
		c += 2;
	}
	
	if (reads.size()) {
		string vname;
		make_volume_name(wrk_dir, vid++, vname);
		reads.dump_pac(vname.c_str());
	}
	
	dump_reads_info(wrk_dir, vid, num_reads);
}
