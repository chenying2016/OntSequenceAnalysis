#ifndef PDB_AUX_H
#define PDB_AUX_H

#include <string>

#include "packed_db.h"

#define DEFAULT_VOLUME_SIZE 2000000000

void
make_volume_name(const char* wrk_dir, const int vid, std::string& name);

//void
//dump_num_volumes(const char* wrk_dir, const int num_vols);

void
dump_reads_info(const char* wrk_dir, const int num_vols, const int num_reads);

int
load_num_volumes(const char* wrk_dir);

int 
load_num_reads(const char* wrk_dir);

void 
validate_platform(const int p);

void
pack_reads(int argc, char* argv[], PackedDB& reads);

void
load_packed_reads(const char* wrk_dir, PackedDB& reads);

#endif // PDB_AUX_H
