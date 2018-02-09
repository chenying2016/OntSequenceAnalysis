#include "../klib/kseq.h"
#include "../klib/kstring.h"
#include "../common/ontcns_aux.h"

#include <stdio.h>

KSEQ_DECLARE(gzFile)

static void
split_one_read(FILE* out, kstring_t* read, int* read_id, size_t min_length, size_t max_length)
{
	size_t from = 0, to = 0;
	size_t read_size = kstr_size(*read);
	while (from < read_size) {
		to = from + max_length;
		to = OC_MIN(to, read_size);
		if (to - from >= min_length) {
			fprintf(out, ">%d\n", *read_id);
			++(*read_id);
			for (size_t i = from; i < to; ++i) {
				fprintf(out, "%c", kstr_A(*read, i));
			}
			fprintf(out, "\n");
		}
		from = to;
	}
}

int main(int argc, char* argv[])
{
	if (argc != 5) {
		fprintf(stderr, "USAGE:\n");
		fprintf(stderr, "%s min-length max-length input output\n", argv[0]);
		return 1;
	}
	size_t min_length = atoll(argv[1]);
	size_t max_length = atoll(argv[2]);
	const char* input = argv[3];
	const char* output = argv[4];
	int id = 0;
	DFOPEN(out, output, "w");
	
	DGZ_OPEN(in, input, "r");
    kseq_t* read = kseq_init(in);
    while (kseq_read(read) >= 0) {
		split_one_read(out, &read->seq, &id, min_length, max_length);
    }
    GZ_CLOSE(in);
    kseq_destroy(read);
	
	FCLOSE(out);
}
