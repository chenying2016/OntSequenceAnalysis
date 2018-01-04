#!/bin/bash

### parse parameters
while read LINE ; do
    eval "${LINE}"
done < $1

echo ${PROJECT}
echo ${THREADS}
echo ${PB_READ_LIST}
echo ${ONT_READ_LIST}
echo ${OVLP_FAST_OPTIONS}
echo ${OVLP_SENSITIVE_OPTIONS}
echo ${CNS_FAST_OPTIONS}
echo ${CNS_SENSITIVE_OPTIONS}
echo ${NUM_ITER}

if [ -n "${PB_READ_LIST}" ]; then
    READ_LIST="-p ${PB_READ_LIST}"
    if [ -n "${ONT_READ_LIST}" ]; then
        READ_LIST="${READ_LIST} -n ${ONT_READ_LIST}"
    fi
else
    if [ ! -n "${ONT_READ_LIST}" ]; then
        echo "Read List Is Empty!"
        exit 1;
    else
        READ_LIST="-n ${ONT_READ_LIST}"
    fi
fi

if [ ! -d ${PROJECT} ]; then
    mkdir -p ${PROJECT}
fi
cd ${PROJECT}
PROJECT_DIR=$(pwd)

for ((i=1;i<=${NUM_ITER};i++))
do
    DIR_NAME="WorkDirIter${i}"
    WRK_DIR="${PROJECT_DIR}/${DIR_NAME}"
    CNS_READS="${WRK_DIR}/${DIR_NAME}Cns.fasta"
    UN_CNS_READS="${WRK_DIR}/${DIR_NAME}UnCns.fasta"


    if [ ! -d ${WRK_DIR} ]; then
        mkdir -p ${WRK_DIR}
    fi

    if [ $i -eq 1 ]; then
        ITER_CAN_OPTIONS="-t ${THREADS} ${OVLP_SENSITIVE_OPTIONS}"
        ITER_CNS_OPTIONS="-t ${THREADS} ${CNS_SENSITIVE_OPTIONS}"
    else
        ITER_CAN_OPTIONS="-t ${THREADS} ${OVLP_FAST_OPTIONS}"
        ITER_CNS_OPTIONS="-t ${THREADS} ${CNS_FAST_OPTIONS}"
    fi

    if [ $i -eq ${NUM_ITER} ]; then
        ITER_CNS_OPTIONS="${ITER_CNS_OPTIONS} -f 0"
    else
        ITER_CNS_OPTIONS="${ITER_CNS_OPTIONS} -f 1"
    fi
    export ONTCNS_CAN_OPTIONS=${ITER_CAN_OPTIONS}
    export ONTCNS_CNS_OPTIONS=${ITER_CNS_OPTIONS}

    CMD="OntCns2Cns.sh -w ${WRK_DIR} -r ${CNS_READS} -u ${UN_CNS_READS} ${READ_LIST}"
    OntCns2RunCmd.sh ${CMD}
    if [ $? -ne 0 ]; then
        exit 1;
    fi

    if [ $? -ne ${NUM_ITER} ]; then
        NEXT_READ_LIST="${WRK_DIR}/ReadList.txt"
        echo "${CNS_READS}" > "${NEXT_READ_LIST}"
        echo "${UN_CNS_READS}" >> "${NEXT_READ_LIST}"
        READ_LIST="-p ${NEXT_READ_LIST}"
    fi
done
cd ..
