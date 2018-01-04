#include "argument.h"
#include "pod_darr.h"
#include "smart_assert.h"

#include <fstream>

using namespace std;

void
validate_arguments(std::vector<Argument*>& args)
{
    std::size_t n = args.size(), i;
    for (i = 0; i < n; ++i) args[i]->validate();
}

void 
parse_arguments(std::vector<Argument*>& args, int argc, char* argv[])
{
    int i = 0, c;
    std::size_t n = args.size(), k;
    std::string name;
    while (argc) {
        if (argv[i][0] != '-') {
            mc_error << "an option name must be start with '-': " << argv[i] << eolog;
        }
        name = argv[i] + 1;
        for (k = 0; k != n; ++k) {
            if (name == args[k]->arg_name()) break;
        }
        if (k == n) mc_error << "unrecognised option '" << argv[i] << "'" << eolog;
        ++i;
        --argc;
        c = args[k]->parse_argument(argc, argv + i);
        i += c;
        argc -= c;
    }
}

void
print_prog_usage(const char* prog, std::vector<Argument*>& args, std::ostream& out)
{
	out << "USAGE:" << "\n"
		<< prog << " [options]" << "\n"
		<< "\n"
		<< "options:" << "\n";
	std::size_t n = args.size(), i;
	for (i = 0; i < n; ++i) args[i]->print_usage(out);
}

void 
print_arguments(std::vector<Argument*>& args, std::ostream& out)
{
	std::size_t n = args.size(), i;
	for (i = 0; i < n; ++i) args[i]->print_value(out);
}

void
load_argv_from_config_file(const char* config_file_path, PODArray<char>& cmd, int& argc, PODArray<char*>& argv)
{
	string line;
	PODArray<int> argv_idx;
	dsopen(ifstream, in, config_file_path, ios::in);
	while (getline(in, line)) {
		if (line[0] == '#') continue;
		if (line.size() == 0) continue;
		r_assert(line[0] == '-')(line);
		argv_idx.push_back(cmd.size());
		cmd.push_back(line.c_str(), line.size());
		cmd.push_back(' ');
	}
	sclose(in);
	
	argv.resize(argv_idx.size());
	char* p = cmd.data();
	for (idx i = 0; i < argv_idx.size(); ++i) {
		argv[i] = p + argv_idx[i];
	}
	argc = argv_idx.size();
}

void 
parse_arguments_from_config_file(const char* config_file_path, std::vector<Argument*>& args)
{
	PODArray<char> cmd;
	int argc;
	PODArray<char*> argv;
	load_argv_from_config_file(config_file_path, cmd, argc, argv);
	parse_arguments(args, argc, argv.data());
}
