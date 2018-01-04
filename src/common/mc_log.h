#ifndef MC_LOG_H
#define MC_LOG_H

#include "defs.h"

#include <iostream>
#include <sstream>
#include <string>

#include <cstdlib>
#include <ctime>
#include <cstring>

enum ELogRank
{
    eINFO = 0,
    eWARN = 1,
    eERROR = 2
};

struct MCEpilog
{
	const char* file;
    const char* func;
    int line;

    MCEpilog(const char* _file, const char* _func, int l) : file(_file), func(_func), line(l) {}
};

class MCLogger
{
public:
    MCLogger(const ELogRank& r) : level(r) {}

    std::ostream& stream() {
        return out[level];
    }

    static void get_current_time(char now[]) {
        time_t ltime;
        time(&ltime);
        char* ltime_str = ctime(&ltime);
        std::size_t n = strlen(ltime_str) - 1;
        strncpy(now, ltime_str, n);
        now[n] = '\0';
    }

    void prolog(std::ostream& os) {
        char now[256];
        get_current_time(now);
        os << "[" << now << "] ";
        switch (level) {
            case eINFO:
                os << "INFO";
                break;
            case eWARN:
                os << "WARN";
                break;
            case eERROR:
                os << "ERROR";
                break;
        }
        os << ": ";
    }

    void epilog(const MCEpilog& e) {
        std::ostream& os = std::cout;
        prolog(os);
        os << out[level].str();
        out[level].str("");
        if (e.func) {
            os << " (" << e.file << ", " << e.func << ", " << e.line << ")";
        }
        os << "\n";
        if (level == eERROR) exit(1);
    }

private:
    ELogRank level;
    std::ostringstream out[3];
};

extern MCLogger mc_log;
extern MCLogger mc_warn;
extern MCLogger mc_error;
#define eolog MCEpilog(__FILE__, __func__, __LINE__)
#define eoelog MCEpilog(NULL, NULL, 0)

#define decl_output(type) \
    MCLogger& operator << (MCLogger& log, const type & t);

decl_output(char)
decl_output(char*)
decl_output(std::string)
decl_output(i32)
decl_output(u32)
decl_output(i64)
decl_output(u64)
decl_output(double)
decl_output(MCEpilog)

#endif // MC_LOG_H
