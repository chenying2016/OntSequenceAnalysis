#ifndef RM_WORKER_H
#define RM_WORKER_H

#include "../common/map_options.h"

void
rm_main(MapOptions* options, const char* reads_path, const char* reference_path, const char* output);

#endif // RM_WORKER_H
