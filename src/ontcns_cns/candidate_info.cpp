#include "candidate_info.h"

std::ostream& operator << (std::ostream& out, const EGappedCandidate& can)
{
    out << "q: "
        << "[" << can.qid << ", " << can.qdir << ", " << can.qext << ", " << can.qsize << "] "
        << "t: "
        << "[" << can.sid << ", " << can.sdir << ", " << can.sext << ", " << can.ssize << "] "
        << "score = " << can.score << "\n";
	return out;
}
