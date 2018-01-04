#include "pdb_aux.h"

#include <fstream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "aux_tools.h"
#include "mc_log.h"
#include "fasta_reader.h"

using namespace std;

void
validate_wrk_dir_path(string& wrk_dir)
{
	if (wrk_dir.empty()) {
		mc_error << "empty working folder." << eolog;
	}
	if (wrk_dir[wrk_dir.size() - 1] != '/') wrk_dir += '/';
}

void
validate_wrk_dir_exist(const string& wrk_dir)
{
	if (access(wrk_dir.c_str(), F_OK) == -1) {
		if (mkdir(wrk_dir.c_str(), 0755) == -1) {
			mc_error << "failed to create folder " << wrk_dir << eolog;
		}
	}
}

void
make_volume_name(const char* wrk_dir, const int vid, string& name)
{
	name = wrk_dir;
	validate_wrk_dir_path(name);
	
	std::ostringstream os;
	os << "vol" << vid;
	name += os.str();
}

void
make_nv_name(const char* wrk_dir, string& name)
{
	name = wrk_dir;
	validate_wrk_dir_path(name);
	name += "reads_info.txt";
}

void
dump_reads_info(const char* wrk_dir, const int num_vols, const int num_reads)
{
	string file_name;
	make_nv_name(wrk_dir, file_name);
	dsopen(ofstream, out, file_name.c_str(), ios::out);
	out << num_vols << "\n"
		<< num_reads << "\n";
	sclose(out);
}

int
load_num_volumes(const char* wrk_dir)
{
	string file_name;
	make_nv_name(wrk_dir, file_name);
	dsopen(ifstream, in, file_name.c_str(), ios::in);
	int num_vols;
	in >> num_vols;
	sclose(in);
	return num_vols;
}

int 
load_num_reads(const char* wrk_dir)
{
	string file_name;
	make_nv_name(wrk_dir, file_name);
	dsopen(ifstream, in, file_name.c_str(), ios::in);
	int num_vols, num_reads;
	in >> num_vols >> num_reads;
	sclose(in);
	return num_reads;
}

void 
validate_platform(const int p)
{
	if (p != TECH_PACBIO && p != TECH_NANOPORE) {
		mc_error << "invalid platform option: " << p << eolog;
	}
}

static void
pack_one_file(const char* path,
			  const int platform,
			  PackedDB& reads)
{
	FastaReader reader(path);
	Sequence read;
	while (1) {
		idx s = reader.read_one_seq(read);
		if (s == -1) break;
		reads.add_one_raw_seq(read, platform);
	}
}

static void
pack_one_list(const char* list_path,
			  const int platform,
			  PackedDB& reads)
{
	dsopen(ifstream, in, list_path, ios::in);
	string line;
	while (getline(in, line)) {
		pack_one_file(line.c_str(), platform, reads);
	}
	sclose(in);
}

void
pack_reads(int argc, char* argv[], PackedDB& reads)
{
	int c = 0;
	while (c < argc) {
		const char* list_path = argv[c];
		int platform = atoi(argv[c+1]);
		validate_platform(platform);
		pack_one_list(list_path, platform, reads);
	}
}

void
load_packed_reads(const char* wrk_dir, PackedDB& reads)
{
	int num_vols = load_num_volumes(wrk_dir);
	for (int i = 0; i < num_vols; ++i) {
		string vname;
		make_volume_name(wrk_dir, i, vname);
		PackedDB pdb;
		pdb.load(vname.c_str());
		reads.merge(pdb);
	}
}
