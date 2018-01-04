#ifndef INDEL_DETECTION_PACKED_DB_H
#define INDEL_DETECTION_PACKED_DB_H

#include <vector>
#include <cstring>
#include <string>

#include "defs.h"
#include "pod_darr.h"
#include "sequence.h"
#include "smart_assert.h"

class PackedDB
{
public:
    struct SeqIndex
    {
        idx offset;
        idx size;
		idx header_offset;
		int tech;
        SeqIndex(idx o = 0, idx s = 0, idx h = 0, int p = TECH_PACBIO)
			: offset(o), size(s), header_offset(h), tech(p) {}
    };
	
	idx size2bytes(const idx s) {
		return (s + 3) >> 2;
	}

public:
    static void set_pac_char(u8* p, const idx i, const u8 c)
    {
        p[i >> 2] |= c << ((~i&3)<<1);
    }
    static u8 get_pac_char(const u8* p, const idx i)
    {
        u8 c = p[i >> 2] >> ((~i&3)<<1)&3;
        return c;
    }
    void generate_pac_name(const char* prefix, std::string& pac_name)
    {
        pac_name = prefix;
        pac_name += ".pac";
    }
    void generate_idx_name(const char* prefix, std::string& idx_name)
    {
        idx_name = prefix;
        idx_name += ".idx";
    }

public:
    PackedDB() : m_pac(NULL), m_db_size(0), m_max_db_size(0) {}
    ~PackedDB() { destroy(); }
	
	void merge(PackedDB& rpdb);

    void set_char(const idx i, const u8 c)
    {
#ifndef NDEBUG
        r_assert(i < m_max_db_size)(i)(m_db_size);
#endif
        set_pac_char(m_pac, i, c);
    }
    u8 get_char(const idx i) const
    {
#ifndef NDEBUG
        r_assert(i < m_db_size)(i)(m_db_size);
#endif
        return get_pac_char(m_pac, i);
    }
	
	int seq_tech(const idx i) const {
		return m_seq_idx[i].tech;
	}

    void destroy()
    {
        if (m_pac) delete[] m_pac;
        m_pac = NULL;
        m_db_size = 0;
        m_max_db_size = 0;
        m_seq_idx.clear();
    }

    void clear()
    {
        if (m_db_size)
        {
            idx pac_bytes =size2bytes(m_db_size);
            memset(m_pac, 0, pac_bytes);
        }
        m_db_size = 0;
        m_seq_idx.clear();
		m_headers.clear();
    }

    void alloc(const idx size)
    {
        destroy();
        m_db_size = 0;
        m_max_db_size = size;
        const u64 bytes = size2bytes(size);
        m_pac = new u8[bytes];
        memset(m_pac, 0, bytes);
    }

    idx size() const
    {
        return m_db_size;
    }

    idx num_seqs() const
    {
        return m_seq_idx.size();
    }

    idx seq_offset(const idx i) const
    {
#ifndef NDEBUG
        const idx ns = num_seqs();
        r_assert(i < ns)(i)(ns);
#endif
        return m_seq_idx[i].offset;
    }

    idx seq_size(const idx i) const
    {
#ifndef NDEBUG
        const idx ns = num_seqs();
        r_assert(i < ns)(i)(ns);
#endif
        return m_seq_idx[i].size;
    }
	
	const char* seq_header(const idx i) const
	{
		return m_headers.data() + m_seq_idx[i].header_offset;
	}

    idx offset_to_id(const idx offset) const
    {
#ifndef NDEBUG
        r_assert(offset < m_db_size)(offset)(m_db_size);
#endif
        idx ns = num_seqs();
        idx left = 0, mid = 0, right = ns;
        while(left < right)
        {
            mid = (left + right) >> 1;
            if (offset >= m_seq_idx[mid].offset)
            {
                if (mid == ns - 1) break;
                if (offset < m_seq_idx[mid + 1].offset) break;
                left = mid + 1;
            }
            else
            {
                right = mid;
            }
        }
        return mid;
    }

	void add_one_encode_seq(const char* seq, const idx size);
    void add_one_raw_seq(const char* seq, const idx size);
	void add_one_raw_seq(Sequence& sequence, int platform);
    void get_sequence(const idx from, const idx to, const bool fwd, char* seq) const;
    void get_sequence(const idx rid, const bool fwd, char* seq) const;
    void get_sequence(const idx rid, const idx from, const idx to, const bool fwd, char* seq) const;
    void get_raw_sequence(const idx rid, const idx from, const idx to, const bool fwd, char* seq) const;
    void dump_pac(const char* path);
    void load_pac(const char* path);
	void load_index(const char* path);
    void load_fasta(const char* path);
	
	void load(const char* path);

private:
    void check_and_realloc(const idx added_size);

private:
    u8*                   m_pac;
    idx                   m_db_size;
    idx                   m_max_db_size;
    std::vector<SeqIndex>   m_seq_idx;
	PODArray<char>			m_headers;
};

enum DBType
{
    eFasta,
    eFastq,
	ePac,
    eEmptyFile,
    eUnknown
};

DBType
detect_db_type(const char* path);

#endif //INDEL_DETECTION_PACKED_DB_H
