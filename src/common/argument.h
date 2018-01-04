#ifndef ARGUMENT_H
#define ARGUMENT_H

#include "aux_tools.h"

enum EArgumentType
{
    eInteger,
    eReal,
    eString
};

class Argument
{
public:
    Argument(const char* n, const char* d, bool optnl)
        : name(n),
          desc(d),
          optional(optnl) {
            is_set = false;
          }

    void enable_set() {
        is_set = true;
    }

    void validate_set() {
        if ((!optional) && (!is_set)) {
            mc_error << "argument to option '" << name << "' must be specified." << eolog;
        }
    }

    const std::string& arg_name() const {
        return name;
    }

    const std::string& arg_desc() const {
        return desc;
    }

    virtual ~Argument() {}

    virtual void validate() = 0;
    
    virtual int parse_argument(int argc, char* argv[]) = 0;

    virtual void print_usage(std::ostream& out) = 0;
	
	virtual void print_value(std::ostream& out) const = 0;

protected:
    const std::string name;
    const std::string desc;
    bool optional;
    bool is_set;
};

template <typename T>
class NumericArgument : public Argument
{
public:
    NumericArgument(const char* n, 
            const char* d, 
            bool optnl, 
            EArgumentType type,
            bool has_dv = false, 
            T _dv = 0, 
            bool has_min = false, 
            T _minv = 0, 
            bool has_max = false, 
            T _maxv = 0)
        : Argument(n, d, optnl),
          has_default(has_dv),
          dv(_dv),
          has_minimum(has_min),
          minv(_minv),
          has_maximum(has_max),
          maxv(_maxv),
          arg_type(type) {
              if (has_default) v = dv;
          }

    virtual ~NumericArgument() {}

    virtual void validate() {
        validate_set();
        if (has_minimum && has_maximum) {
            if (v < minv || v > maxv) {
                mc_error << "argument to option '" << arg_name() << "' must be in range [" << minv << " .. " << maxv << "]: " << v << eolog;
            }
        } else if (has_minimum) {
            if (v < minv) {
                mc_error << "argument to option '" << arg_name() << "' must be >= " << minv << ": " << v << eolog;
            }
        } else if (has_maximum) {
            if (v > maxv) {
                mc_error << "argument to option '" << arg_name() << "' must be <= " << maxv << ": " << v << eolog;
            }
        }
    }

    virtual void print_usage(std::ostream& out) {
        out << "  " << "-" << arg_name() << "<";
        switch (arg_type) {
            case eInteger:
                out << "Integer";
                break;
            case eReal:
                out << "Real";
                break;
			default:
				break;
        }
        if (has_minimum && has_maximum) {
            out << ", [" << minv << " .. " << maxv << "]";
        } else if (has_minimum) {
            out << ", >= " << minv;
        } else if (has_maximum) {
            out << ", <= " << maxv;
        }
        out << ">";
        out << "\n";
        out << "    " << arg_desc() << "\n";
        if (has_default) out << "    " << "default: " << dv << "\n";
    }

    virtual int parse_argument(int argc, char* argv[]) {
        if (!argc) {
            mc_error << "argument to option '" << arg_name() << "' is missing" << eolog;
        }
        std::istringstream in(argv[0]);
        in >> v;
        if (!in) mc_error << "fail to parse argument '" << argv[0] << "' to option '" << arg_name() << "'" << eolog;
        enable_set();
        return 1;
    }

    T value() const {
        return v;
    }
	
	virtual void print_value(std::ostream& out) const {
		out << arg_name() << ":\t\t" << value() << "\n";
	}

private:
    T v;
    bool has_default;
    T dv;
    bool has_minimum;
    T minv;
    bool has_maximum;
    T maxv;
    EArgumentType arg_type;
};

class StringArgument : public Argument
{
public:
    StringArgument(const char* n, const char* d, bool optnl, bool has_dv = false, const char* _dv = NULL)
        : Argument(n, d, optnl),
          has_default(has_dv),
          dv(_dv) {
              if (has_default) v = dv;
          }

    virtual ~StringArgument() {}

    virtual void validate() {
        validate_set();
    }

    virtual void print_usage(std::ostream& out) {
        out << "  " << "-" << arg_name() << "<string>" << "\n";
        out << "    " << arg_desc() << "\n";
        if (has_default) {
            out << "    " << "default: " << dv << "\n";
        }
    }

    virtual int parse_argument(int argc, char* argv[]) {
        if (!argc) mc_error << "argument to option '" << arg_name() << "' is missing" << eolog;
        v = argv[0];
		enable_set();
        return 1;
    }

    const char* value() const {
        return v;
    }
	
	virtual void print_value(std::ostream& out) const {
		out << arg_name() << ":\t\t";
		if (v) out << v;
		out << "\n";
	}

private:
    const char* v;
    bool has_default;
    const char* dv;
};

class FlagArgument : public Argument
{
public:
    FlagArgument(const char* n, const char* d, bool has_dv, bool _dv) 
        : Argument(n, d, true),
          has_default(has_dv),
          dv(_dv) {
              if (has_default) flag_is_set = dv;
          }

    virtual ~FlagArgument() {}

    virtual void validate() {}

    virtual int parse_argument(int, char*[]) {
        flag_is_set = true;
		enable_set();
        return 0;
    }

    virtual void print_usage(std::ostream& out) {
        out << "  " << "-" << arg_name() << "\n"
            << "    " << arg_desc() << "\n";
    }
	
	bool value() const {
		return flag_is_set;
	}
	
	virtual void print_value(std::ostream& out) const {
		out << arg_name() << ":\t\t";
		if (value()) out << "enabled";
		else out << "disable";
		out << "\n";
	}
	
private:
    bool flag_is_set;
    bool has_default;
    bool dv;
};

typedef NumericArgument<i32> Int4Argument;
typedef NumericArgument<i64> Int8Argument;
typedef NumericArgument<double> RealArgument;

#include <vector>

void
validate_arguments(std::vector<Argument*>& args);

void 
parse_arguments(std::vector<Argument*>& args, int argc, char* argv[]);

void
print_prog_usage(const char* prog, std::vector<Argument*>& args, std::ostream& out);

void 
print_arguments(std::vector<Argument*>& args, std::ostream& out);

void 
parse_arguments_from_config_file(const char* config_file_path, std::vector<Argument*>& args);

#endif // ARGUMENT_H
