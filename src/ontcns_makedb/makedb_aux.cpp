#include "makedb_aux.h"

#include <fstream>

#include "../common/aux_tools.h"
#include "../common/mc_log.h"

using namespace std;

void
make_volume_name(const char* wrk_dir, const int vid, std::string& name)
{
	name = wrk_dir;
	if (name[ name.size() - 1] != '/') name += '/';
	
	std::ostringstream os;
	os << "vol" << vid;
	name += os.str();
}

void
make_nv_name(const char* wrk_dir, string& name)
{
	name = wrk_dir;
	if (name.empty()) {
		mc_error << "empty working folder" << eolog;
	}
	if (name[name.size() - 1] != '/') name += '/';
	name += "num_vols.txt";
}

void
dump_num_vols(const char* wrk_dir, const int num_vols)
{
	string file_name;
	make_nv_name(wrk_dir, file_name);
	dsopen(ofstream, out, file_name.c_str(), ios::out);
	out << num_vols << "\n";
	sclose(out);
}

int
load_num_vols(const char* wrk_dir)
{
	string file_name;
	make_nv_name(wrk_dir, file_name);
	dsopen(ifstream, in, file_name.c_str(), ios::in);
	int num_vols;
	in >> num_vols;
	sclose(in);
	return num_vols;
}
