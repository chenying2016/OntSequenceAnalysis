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
	./common/kstring.c \
	./common/makedb_aux.c \
	./common/map_aux.c \
	./common/map_options.c \
	./common/nst_nt4_table.c \
	./common/packed_db.c \
	./common/sequence.c \
	./common/results_writer.c \
	./lookup_table/hash_list_bucket_sort.c \
	./lookup_table/lookup_table.c \
	./word_finder/chain_dp.c \
	./word_finder/word_finder_aux.c \
	./word_finder/word_finder.c

SRC_INCDIRS  := common \

SUBMAKEFILES := ./test/main.mk \
	./makedb/main.mk \
	./pm_one_volume/main.mk \
	./pairwise_mapping/main.mk
