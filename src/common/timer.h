#ifndef TIMER_H
#define TIMER_H

#include "mc_log.h"
#include <sys/time.h>

struct Timer
{
    struct timeval start;
    struct timeval end;

    void go() { gettimeofday(&start, NULL); }
    void stop() { gettimeofday(&end, NULL); }

    double elapsed() {
        return end.tv_sec - start.tv_sec + 1.0 * (end.tv_usec - start.tv_usec) / 1000000;
    }

    void elapsed(char dur[]) {
        double d = elapsed();
        sprintf(dur, "%.2lf", d);
    }
};

struct DynamicTimer
{
	DynamicTimer(const char* func) : m_func(func) {
		mc_log << "'" << m_func << "' begins" << eoelog;
		timer.go();
	}
	
	~DynamicTimer() {
		timer.stop();
		char dur[256];
		timer.elapsed(dur);
		mc_log << "'" << m_func << "' takes " << dur << " secs." << eoelog;
	}
	
private:
	const char* m_func;
	Timer timer;
};

#endif // TIMER_H


