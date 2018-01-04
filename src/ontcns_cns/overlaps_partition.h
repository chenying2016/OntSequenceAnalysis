#ifndef OVERLAPS_PARTITION_H
#define OVERLAPS_PARTITION_H

#include "cns_options.h"

#include <string>
#include <vector>

void
make_partition_file_name(const char* prefix, const int part, std::string& name);

int
partition_candidates(const CnsOptions& options);

#endif // OVERLAPS_PARTITION_H
