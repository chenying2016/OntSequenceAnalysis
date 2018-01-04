#ifndef RECORD_READER_H
#define RECORD_READER_H

#include "aux_tools.h"
#include "buffer_line_reader.h"
#include <sstream>
#include <string>

class RecordReader: public mc_noncopyable
{
public:
    RecordReader(const char* record_path):
        line_reader_(record_path) {}
	
	template <typename T>
	bool read_one_record(T& v) {
		if (++line_reader_) {
			const BufferLineReader::OneDataLine& line_data = line_reader_.get_line();
			line.assign(line_data.begin(), line_data.end());
			ins.str( line );
			ins.clear();
			v.read_from(ins);
			return 1;
		} else {
			return 0;
		}
	}

private:
	std::string line;
    std::istringstream ins;
    BufferLineReader line_reader_;
};

#endif // RECORD_READER_H
