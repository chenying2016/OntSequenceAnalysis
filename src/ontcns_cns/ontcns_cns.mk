ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := OntCns2Consensus
SOURCES  := main.cpp \
	cns_options.cpp \
	overlap_pool.cpp \
	overlaps_partition.cpp \
	candidate_info.cpp \
	mc_consensus.cpp \
	consensus_aux.cpp \
	cns_stage2.cpp

SRC_INCDIRS  := . 

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lontcns
TGT_PREREQS := libontcns.a

SUBMAKEFILES :=
