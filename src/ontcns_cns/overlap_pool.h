#ifndef OVERLAP_POOL_H
#define OVERLAP_POOL_H

#include "candidate_info.h"
#include "../common/smart_assert.h"

struct EditScriptWorker
{
    static u8 pack_edit_script(const u8 type, const u8 num) {
        r_assert(num <= kMaxNum);
        return type | num;
    }

    static void unpack_edit_script(const u8 es, int& type, int& num) {
        type = es & kTypeMask;
        num = es & kNumMask;
    }

    static const int kMaxNum    = 63;
    static const u8 kOpM        = 64;
    static const u8 kOpI        = 128;
    static const u8 kOpD        = 192;
    static const u8 kTypeMask   = 192;
    static const u8 kNumMask    = 63;
};

void
add_one_edit_segment(const u8 type, int num, std::vector<u8>& scripts);

void
build_edit_scripts(const std::string& qaln, const std::string& taln, std::vector<u8>& scripts);

class OverlapPool
{
public:
    OverlapPool() :
        edit_scripts(new u8[kPoolSize]),
        cur(0) {
            std::fill(next_avail, next_avail + kMaxIdx, -1);
            std::fill(add_cnt, add_cnt + kMaxIdx + 2, 0);
            last_valid_pci_idx = 0;
        }

    ~OverlapPool() {
        //for (int i = 0; i < kMaxIdx + 2; ++i) {
        //    std::cout << "add[ " << i << "] = " << add_cnt[i] << "\n";
        //}
        //mc_log << "overlap pool used size: " << cur << eolog;
        destroy();
    }

    void destroy() {
        delete[] edit_scripts;
        edit_scripts = NULL;
        cur = 0;
    }

    idx used_size() const {
        return cur;
    }

    void add_scripts(std::vector<u8>& scripts, int& es_start);

    void get_scripts(EGappedCandidate& c, std::vector<u8>& scripts);

    const u8* get_aln_scripts(const idx offset) const {
        return edit_scripts + offset;
    }

    void dealloc_invalid_scripts(int end_ci_id, int max_read_id, EGCInfo* p);

private:
    static const idx kPoolSize      = (idx)1 * 1024 * 1024 * 1024;
    static const idx kMaxScriptSize = 20000;
    static const idx kInc           = 1000;
    static const idx kMaxIdx        = 20;

    idx next_avail[kMaxIdx];
    int add_cnt[kMaxIdx + 2];
    u8* edit_scripts;
    idx cur;
    idx last_valid_pci_idx;
};

#endif // OVERLAP_POOL_H
