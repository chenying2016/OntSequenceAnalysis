#!/bin/bash

### parse arguments
while read LINE ; do
    eval "${LINE}"
done < $1

echo ${PROJECT}
echo ${THREADS}
echo ${ONT_READ_LIST}
echo ${OVLP_FAST_OPTIONS}
echo ${OVLP_SENSITIVE_OPTIONS}
echo ${CNS_FAST_OPTIONS}
echo ${CNS_SENSITIVE_OPTIONS}
echo ${NUM_ITER}

READ_LIST="-n ${ONT_READ_LIST}"

if [ ! -d ${PROJECT} ]; then
    mkdir -p ${PROJECT}
fi
cd ${PROJECT}
PROJECT_DIR=$(pwd)

for ((i=1;i<=${NUM_ITER};i++))
do
    DIR_NAME="cns_iter${i}"
    WRK_DIR="${PROJECT_DIR}/${DIR_NAME}"
    CNS_READS="${WRK_DIR}/cns.fasta"
    UNCNS_READS="${WRK_DIR}/raw.fasta"

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

    CMD="oc2cns.sh -w ${WRK_DIR} -r ${CNS_READS} -u ${UNCNS_READS} ${READ_LIST}"
    oc2cmd.sh ${CMD}
    if [ $? -ne 0 ]; then
        exit 1;
    fi

    if [ $i -ne ${NUM_ITER} ]; then
        NEXT_READ_LIST="${WRK_DIR}/ReadList.txt"
        echo "${CNS_READS}" > "${NEXT_READ_LIST}"
        echo "${UNCNS_READS}" >> "${NEXT_READ_LIST}"
        READ_LIST="-n ${NEXT_READ_LIST}"
    fi
done
cd ..
