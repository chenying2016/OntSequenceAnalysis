#ifndef MAKEDB_AUX_H
#define MAKEDB_AUX_H

#include <string>

void
make_volume_name(const char* wrk_dir, const int vid, std::string& name);

void
dump_num_vols(const char* wrk_dir, const int num_vols);

int
load_num_vols(const char* wrk_dir);

#endif // MAKEDB_AUX_H
