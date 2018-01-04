#include <iostream>
#include "mc_consensus.h"

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
	CnsOptions options;
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
	
	cns_main(options);
	return 0;
}
