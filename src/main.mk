ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libontcns.a

SOURCES      := \
		common/argument.cpp \
		common/buffer_line_reader.cpp \
		common/cns_seq.cpp \
		common/defs.cpp \
		common/fasta_reader.cpp \
		common/gapped_candidate.cpp \
		common/lookup_table.cpp \
		common/lookup_table_bucket_sort.cpp \
		common/m4record.cpp \
		common/mc_log.cpp \
		common/mc_ostream.cpp \
		common/output_stream.cpp \
		common/packed_db.cpp \
		common/pdb_aux.cpp \
		common/sequence.cpp \
		tasc/align_tags.cpp \
		tasc/cbcns.cpp \
		tasc/cns_aux.cpp \
		tasc/soa.cpp \
		word_finder/chain_dp.cpp \
		word_finder/word_finder_aux.cpp \
		word_finder/word_finder.cpp \
		gapped_align/diff_ex.cpp \
		gapped_align/edlib_ex.cpp \
		gapped_align/ontcns_aligner.cpp \
		gapped_align/results_pool.cpp \
		gapped_align/xdrop_sw.cpp

SRC_INCDIRS  := common \
				word_finder

SUBMAKEFILES := ./ontcns_cns/ontcns_cns.mk \
	./sequence_stats/ontcns_sequence_stats.mk \
	./split_long_reads/ontcns_split_long_reads.mk \
	./ontcns_makedb/ontcns_makedb.mk \
	./ontcns_candidate_detector/OntCns2CanFinder.mk
