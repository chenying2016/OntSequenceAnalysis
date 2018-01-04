#include "options.h"
#include "can_finder.h"

#include <iostream>

#include <cstring>

using namespace std;

int main(int argc, char* argv[])
{
    MapOptions options;
    int r = parse_arguments(argc, argv, options);
    if (r == -1) {
        print_usage();
        return 1;
    }
    if (r == 1) {
        print_usage();
        return 0;
    }
	
	print_options(options);
	
	candidate_main(options);
}
