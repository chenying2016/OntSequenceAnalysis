#include "overlap_pool.h"

using namespace std;

const int EditScriptWorker::kMaxNum;

void add_one_edit_segment(const u8 type, int num, vector<u8>& scripts)
{
    const int MaxSize = EditScriptWorker::kMaxNum;
    int left = num;
    while (left) {
        int n = min(left, MaxSize);
        u8 es = EditScriptWorker::pack_edit_script(type, n);
        scripts.push_back(es);
        left -= n;
    }
}

void
build_edit_scripts(const string& qaln, const string& taln, vector<u8>& scripts)
{
    scripts.clear();
    r_assert(qaln.size() == taln.size());
    int i = 0, j = 0, n = qaln.size();
    while (i < n) {
        if (qaln[i] == '-') {
            j = i + 1;
            while (j < n && qaln[j] == '-') ++j;
            add_one_edit_segment(EditScriptWorker::kOpD, j - i, scripts);
        } else if (taln[i] == '-') {
            j = i + 1;
            while (j < n && taln[j] == '-') ++j;
            add_one_edit_segment(EditScriptWorker::kOpI, j - i, scripts);
        } else {
            j = i + 1;
            while (j < n && qaln[j] != '-' && taln[j] != '-') ++j;
            add_one_edit_segment(EditScriptWorker::kOpM, j - i, scripts);
        }
        i = j;
    }
}

inline bool
get_next_idx(const u8* ed, int& next_idx)
{
    if (!*ed) return 0;
    memcpy( (void*)&next_idx, (void*)(ed + 1), sizeof(int) );
    return 1;
}

inline void
set_next_idx(u8* ed, int next_idx)
{
    if (next_idx == -1) {
        *ed = 0;
        return;
    }
    *ed = 1;
    memcpy( (void*)(ed + 1), (void*)&next_idx, sizeof(int) );
}

inline int
adjust_es_size(int es_size, int inc)
{
    return es_size + inc - (es_size % inc);
}

void
OverlapPool::add_scripts(vector<u8>& scripts, int& es_start)
{
    const int es_size = scripts.size();

    if (es_size >= kMaxScriptSize) {
        es_start = -1;
        add_cnt[kMaxIdx]++;
        return;
    }

    int bid = es_size / kInc;
    ++add_cnt[bid];
    if (next_avail[bid] == -1) {
        if (cur + scripts.size() >= kPoolSize) {
            ++add_cnt[kMaxIdx + 1];
            es_start = -1;
        } else {
            es_start = cur;
            memcpy(edit_scripts + cur, scripts.data(), scripts.size());
            cur += adjust_es_size(es_size, kInc);
        }
    } else {
        es_start = next_avail[bid];
        int next_idx;
        if (get_next_idx(edit_scripts + es_start, next_idx)) {
            next_avail[bid] = next_idx;
        } else {
            next_avail[bid] = -1;
        }
        memcpy(edit_scripts + es_start, scripts.data(), scripts.size());
    }
}

void
OverlapPool::get_scripts(EGappedCandidate& c, vector<u8>& scripts)
{
    const int es_start = c.es_start;
    const int es_size = c.es_size;
    r_assert(es_start % kInc == 0);
    scripts.resize(es_size);
    memcpy(scripts.data(), edit_scripts + es_start, es_size);
    int bid = es_size / kInc;
    set_next_idx(edit_scripts + es_start, next_avail[bid]);
    next_avail[bid] = es_start;
    c.status = EGappedCandidate::kInvalid;
}

void
OverlapPool::dealloc_invalid_scripts(int end_ci_id, int max_read_id, EGCInfo* p)
{
    int num_invalid_ci = 0;
    int num_valid_ci = 0;
    int i;
    for (i = last_valid_pci_idx; i < end_ci_id; ++i) {
        EGCInfo& ci = p[i];
        EGappedCandidate& ec = *ci.p;
        if (ec.qid >= max_read_id || ec.sid >= max_read_id) {
            if (ec.status == EGappedCandidate::kAligned) ++num_valid_ci;
            continue;
        }
        if (ec.status == EGappedCandidate::kAligned) {
            const int bid = ec.es_size / kInc;
            set_next_idx(edit_scripts + ec.es_start, next_avail[bid]);
            next_avail[bid] = ec.es_start;
            ec.status = EGappedCandidate::kValid;
            ++num_invalid_ci;
        }
    }
    int j = end_ci_id - 1;
    for (i = end_ci_id - 1; i >= last_valid_pci_idx; --i) {
        if (p[i].p->status == EGappedCandidate::kAligned) {
            p[j--] = p[i];
        }
    }
    ++j;
    r_assert(end_ci_id - num_invalid_ci == j);
    last_valid_pci_idx = j;

    mc_log << "last_valid_pci_idx: " << last_valid_pci_idx 
           << ", number of valid can_info: " << num_valid_ci
           << ", number of invalid can_info: " << num_invalid_ci
           << eolog;
}
