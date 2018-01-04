#include "mc_log.h"

using namespace std;

MCLogger mc_log(eINFO);
MCLogger mc_warn(eWARN);
MCLogger mc_error(eERROR);

#define def_output(type) \
MCLogger& operator << (MCLogger& log, const type & t) \
{ \
    log.stream() << t; \
    return log; \
}

def_output(char)
def_output(char*)
def_output(string)
def_output(i32)
def_output(u32)
def_output(i64)
def_output(u64)
def_output(double)

MCLogger& operator << (MCLogger& log, const MCEpilog& e)
{
    log.epilog(e);
    return log;
}

