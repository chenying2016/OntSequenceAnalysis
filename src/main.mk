ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libontcns.a

SOURCES      := \
	./common/gapped_candidate.c \
	./common/ontcns_aux.c \
	./common/ontcns_defs.c \
	./common/m4_record.c \
	./common/makedb_aux.c \
	./common/map_aux.c \
	./common/map_options.c \
	./common/nst_nt4_table.c \
	./common/oc_assert.c \
	./common/packed_db.c \
	./common/record_reader.c \
	./common/record_writer.c \
	./common/soa.c \
	./gapped_align/edlib_ex.c \
	./gapped_align/edlib_ex_aux.c \
	./gapped_align/oc_aligner.c \
	./klib/kstring.c \
	./lookup_table/hash_list_bucket_sort.c \
	./lookup_table/lookup_table.c \
	./partition_candidates/pcan_aux.c \
	./tasc/align_tags.c \
	./tasc/cbcns.c \
	./tasc/cns_aux.c \
	./tasc/align_tags.c \
	./word_finder/chain_dp.c \
	./word_finder/word_finder_aux.c \
	./word_finder/word_finder.c

SRC_INCDIRS  := common \

SUBMAKEFILES := ./test/main.mk \
	./makedb/main.mk \
	./pm_one_volume/main.mk \
	./pairwise_mapping/main.mk \
	./partition_candidates/main.mk \
	./consensus/main.mk \
	./reference_mapping/main.mk \
	./assembly/main.mk \
	./sequence_length_stats/main.mk \
	./split_long_reads/main.mk
