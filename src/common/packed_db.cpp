
#include <ios>
#include <fstream>
#include "defs.h"
#include "fasta_reader.h"
#include "packed_db.h"
#include "smart_assert.h"
#include "timer.h"

using namespace std;

static const char* pac_header = "fowfgeufjfwoehugg5926r79";

void 
write_pac_header(std::ostream& out)
{
	int s = strlen(pac_header);
	sbwrite(out, pac_header, s);
}

bool 
check_pac_header(std::istream& in)
{
	int s = strlen(pac_header);
	char buf[256];
	sbread(in, buf, s);
	buf[s] = '\0';
	bool r = (strcmp(buf, pac_header) == 0);
	return r;
}

bool 
check_pac_header(const char* path)
{
	dsopen(ifstream, in, path, ios::binary);
	bool r = check_pac_header(in);
	sclose(in);
	return r;
}

DBType
detect_db_type(const char* path)
{
	if (access(path, F_OK) == -1) {
		mc_error << "file '" << path << "' does not exit" << eolog;
	}
	if (!file_size(path)) {
		return eEmptyFile;
	}
	if (check_pac_header(path)) {
		return ePac;
	}
	
    DBType dbt = eUnknown;
    dsopen(std::ifstream, in, path, std::ios::in);
    std::string line;
    if (!getline(in, line)) dbt = eEmptyFile;
    else if (line[0] == '>') dbt = eFasta;
    else if (line[0] == '@') dbt = eFastq;
    sclose(in);
    return dbt;
}

void
PackedDB::merge(PackedDB& rpdb)
{
	/// index
	const idx ho_start = m_headers.size();
	const idx seq_start = m_db_size;
	for (size_t i = 0; i != rpdb.m_seq_idx.size(); ++i) {
		SeqIndex si = rpdb.m_seq_idx[i];
		si.offset += seq_start;
		si.header_offset += ho_start;
		m_seq_idx.push_back(si);
	}
	
	/// headers
	//m_headers.insert(m_headers.end(), rpdb.m_headers.begin(), rpdb.m_headers.end());
	m_headers.push_back(rpdb.m_headers.begin(), rpdb.m_headers.size());
	
	/// db size
	check_and_realloc(rpdb.m_db_size);
	
	/// pac
	for (idx i = 0; i < rpdb.m_db_size; ++i) {
		u8 c = get_pac_char(rpdb.m_pac, i);
		set_pac_char(m_pac, m_db_size, c);
		++m_db_size;
	}
}

void
PackedDB::check_and_realloc(const idx added_size)
{
    if (m_db_size + added_size > m_max_db_size)
    {
        idx new_size = (m_max_db_size) ? m_max_db_size : 4096;
        while (m_max_db_size + added_size > new_size) new_size *= 2;
        idx bytes = size2bytes(new_size);
        u8* new_pac = new u8[bytes];
        memset(new_pac, 0, bytes);
        memcpy(new_pac, m_pac,(m_db_size + 3) / 4);
        if (m_pac) delete[] m_pac;
        m_pac = new_pac;
        m_max_db_size = new_size;
    }
}

void
PackedDB::add_one_raw_seq(const char *seq, const idx size)
{
    m_seq_idx.push_back(SeqIndex(m_db_size, size));
    check_and_realloc(size);
    const u8* et = get_dna_encode_table();
    for (idx i = 0; i < size; ++i)
    {
        u8 c = seq[i];
        c = et[c];
        if (c > 3) c = rand() & 3;
        set_char(m_db_size, c);
        ++m_db_size;
    }
}

void 
PackedDB::add_one_encode_seq(const char* seq, const idx size)
{
	check_and_realloc(size);
	for (idx i = 0; i < size; ++i)
	{
		u8 c = seq[i];
		set_char(m_db_size, c);
		++m_db_size;
	}
}

void
PackedDB::add_one_raw_seq(Sequence& sequence, int platform)
{
	m_seq_idx.push_back(SeqIndex(m_db_size, sequence.size(), m_headers.size(), platform));
	check_and_realloc(sequence.size());
	const u8* et = get_dna_encode_table();
	const char* seq = sequence.sequence().data();
	for (idx i = 0; i < sequence.size(); ++i)
	{
		u8 c = seq[i];
        c = et[c];
        if (c > 3) c = rand() & 3;
        set_char(m_db_size, c);
        ++m_db_size;
	}
	m_headers.push_back(sequence.header().data(), sequence.header().size());
	m_headers.push_back('\0');
}

void
PackedDB::get_sequence(const idx from, const idx to, const bool fwd, char *seq) const
{
//#ifndef NDEBUG
    r_assert(from <= to)(from)(to);
    r_assert(to <= m_db_size)(to)(m_db_size);
//#endif

    idx pos = 0;
    if (fwd)
         for (idx i = from; i < to; ++i)
         {
             u8 c = get_char(i);
             seq[pos++] = c;
         }
    else
         for (idx i = to - 1; i >= from; --i)
         {
             u8 c = get_char(i);
             c = 3 - c;
             seq[pos++] = c;
         }
    //r_assert(pos <= MAX_SEQ_SIZE)(pos);
}

void
PackedDB::get_sequence(const idx rid, const bool fwd, char *seq) const
{
#ifndef NDEBUG
    const idx ns = num_seqs();
    r_assert(rid < ns)(rid)(ns);
#endif

    const idx s = seq_offset(rid);
    const idx e = s + seq_size(rid);
    get_sequence(s, e, fwd, seq);
}

void
PackedDB::get_sequence(const idx rid, const idx from, const idx to, const bool fwd, char *seq) const
{
#ifndef NDEBUG
    const idx ns = num_seqs();
    r_assert(rid < ns)(rid)(ns);
    const idx ss = seq_size(rid);
    r_assert(from <= to)(from)(to);
    r_assert(to <= ss)(to)(ss);
#endif

    const idx s = seq_offset(rid) + from;
    const idx e = seq_offset(rid) + to;
    get_sequence(s, e, fwd, seq);
}

void
PackedDB::get_raw_sequence(const idx rid, const idx from, const idx to, const bool fwd, char *seq) const
{
    get_sequence(rid, from, to, fwd, seq);
	idx i;
    for (i = 0; i < to - from; ++i)
    {
        u8 c = seq[i];
        c = DNADecoder::decode(c); //"ACGT"[c];
        seq[i] = c;
    }
	seq[i] = '\0';
}

void
PackedDB::dump_pac(const char *path) {
    dsopen(std::ofstream, out, path, std::ios::binary);
	write_pac_header(out);
    const u64 idx_bytes = sizeof(idx);
    // 1) num seqs
    const idx ns = num_seqs();
    sbwrite(out, &ns, idx_bytes);
    // 2) seq idx
    SeqIndex* seq_idx = m_seq_idx.data();
    const u64 idx_size = sizeof(SeqIndex) * ns;
    sbwrite(out, seq_idx, idx_size);
	// 3) header size
	const idx header_size = m_headers.size();
	sbwrite(out, &header_size, idx_bytes);
	// 4 headers
	sbwrite(out, m_headers.data(), header_size);
    // 5) db size
    sbwrite(out, &m_db_size, idx_bytes);
    // 6) pac
    const u64 pac_bytes = size2bytes(m_db_size);
    sbwrite(out, m_pac, pac_bytes);

    sclose(out);
}

void
PackedDB::load_index(const char *path) {
	clear();
    dsopen(std::ifstream, in, path, std::ios::binary);
	if (!check_pac_header(in)) {
		mc_error << "file '" << path << "' is not in pac format" << eolog;
	}
    const u64 idx_bytes = sizeof(idx);
    // 1) num seqs
    idx ns;
    sbread(in, &ns, idx_bytes);
    // 2) seq idx
    m_seq_idx.resize(ns);
    SeqIndex* seq_idx = m_seq_idx.data();
    const u64 idx_size = sizeof(SeqIndex) * ns;
    sbread(in, seq_idx, idx_size);
	// 3) header_size
	idx header_size;
	sbread(in, &header_size, idx_bytes);
	// 4) headers
	m_headers.clear();
	m_headers.resize(header_size);
	sbread(in, m_headers.data(), header_size);
    // 5) db_size
    sbread(in, &m_db_size, idx_bytes);
	sclose(in);
}

void
PackedDB::load_pac(const char *path) {
    clear();
    dsopen(std::ifstream, in, path, std::ios::binary);
	if (!check_pac_header(in)) {
		mc_error << "file '" << " is not in pac format" << eolog;
	}
    const u64 idx_bytes = sizeof(idx);
    // 1) num seqs
    idx ns;
    sbread(in, &ns, idx_bytes);
    // 2) seq idx
    m_seq_idx.resize(ns);
    SeqIndex* seq_idx = m_seq_idx.data();
    const u64 idx_size = sizeof(SeqIndex) * ns;
    sbread(in, seq_idx, idx_size);
	// 3) header_size
	idx header_size;
	sbread(in, &header_size, idx_bytes);
	// 4) headers
	m_headers.clear();
	m_headers.resize(header_size);
	sbread(in, m_headers.data(), header_size);
    // 5) db_size
	idx db_size;
    sbread(in, &db_size, idx_bytes);
	// 6) pac
    check_and_realloc(db_size);
    const u64 pac_bytes = size2bytes(db_size);
    sbread(in, m_pac, pac_bytes);
	m_db_size = db_size;
    sclose(in);
}

void
PackedDB::load_fasta(const char *path)
{
    destroy();
    DBType dbt = detect_db_type(path);
    if (dbt == eEmptyFile) return;
	if (dbt == eUnknown) mc_error << "unknown format file '" << path << "'" << eolog;
    idx db_size = file_size(path);
    if (dbt == eFastq) db_size /= 2;
    db_size += MAX_SEQ_SIZE;
    check_and_realloc(db_size);
    FastaReader reader(path);
    Sequence seq;
	int seq_id = 0;
    while (1)
    {
        idx size = reader.read_one_seq(seq);
        if (size == -1) break;
		check_and_rename_ontcns_header(seq, seq_id++);
		add_one_raw_seq(seq, TECH_PACBIO);
    }
}

void
PackedDB::load(const char* path)
{
	DBType dbt = detect_db_type(path);
	if (dbt == eEmptyFile) return;
	if (dbt == ePac) {
		load_pac(path);
	} else if (dbt == eFasta || dbt == eFastq) {
		load_fasta(path);
	} else {
		mc_error << "unknown format file: '" << path << "'" << eolog;
	}
}
