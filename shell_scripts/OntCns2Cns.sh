#!/bin/bash

PARSE_ARG_SUCCESS=1
WRK_DIR=""
CAN_OPTIONS="${ONTCNS_CAN_OPTIONS}"
CNS_OPTIONS="${ONTCNS_CNS_OPTIONS}"
PB_READS_LIST=""
ONT_READS_LIST=""
READS_LIST=""
CANDIDATES=""
CNS_RESULTS=""
UNCNS_RESULTS=""
PACK_WORKER="OntCns2MakeDB"
CAN_WORKER="OntCns2CanFinder"
CNS_WORKER="OntCns2Consensus"
PACKED_DATA_DIR=""


fPrintUsage()
{
    echo "USAGE:"
    echo $0 -w wrk_dir -r results -u uncorrected_results -p pb_reads_list -n ont_reads_list
}

fValidateOptions()
{
    if [ ! -n "${WRK_DIR}" ]; then
        echo "Working Folder (-w) Is Not Set!"
        PARSE_ARG_SUCCESS=0;
    fi
    if [ ! -n "${READS_LIST}" ]; then
        echo "Read List (-p or -n) Is Empty!"
        PARSE_ARG_SUCCESS=0;
    fi
    if [ ! -n "${CNS_RESULTS}" ]; then
        echo "Corrected Results Path Is Not Specified!"
        PARSE_ARG_SUCCESS=0;
    fi
    if [ ! -n "${UNCNS_RESULTS}" ]; then
        echo "Uncorrected Results Path Is Not Set!"
        PARSE_ARG_SUCCESS=0;
    fi
}

while getopts "w:p:n:r:u:" OPTION; 
do
    case ${OPTION} in
        w)
            WRK_DIR=${OPTARG}
            ;;
        p)
            PB_READS_LIST=${OPTARG}
            ;;
        n)
            ONT_READS_LIST=${OPTARG}
            ;;
        r)
            CNS_RESULTS=${OPTARG}
            ;;
        u)
            UNCNS_RESULTS=${OPTARG}
            ;;
        ?)
            PARSE_ARG_SUCCESS=0
            break;
            ;;
    esac
done

if [ -n "${PB_READS_LIST}" ]; then
    READS_LIST="${PB_READS_LIST} 0"
fi

if [ -n "${ONT_READS_LIST}" ]; then
    READS_LIST="${READS_LIST} ${ONT_READS_LIST} 1"
fi

fValidateOptions;

if [ ${PARSE_ARG_SUCCESS} -eq 0 ]; then
    echo "Parse Argument Failed!"
    fPrintUsage;
    exit 1;
fi

PACK_JOB_FINISHED="${WRK_DIR}/pac.finished"
CAN_JOB_FINISHED="${WRK_DIR}/can.finished"
CNS_JOB_FINISHED="${WRK_DIR}/cns.finished"
if [ -f ${CNS_JOB_FINISHED} ]; then
    echo "Job [${WRK_DIR}] Has Been Finished, Skip It."
    exit 0;
fi

PACKED_DATA_DIR="${WRK_DIR}/PackedData"
echo "Data Dir: ${PACKED_DATA_DIR}"
if [ ! -d ${PACKED_DATA_DIR} ]; then
    mkdir -p ${PACKED_DATA_DIR}
fi
CANDIDATES="${WRK_DIR}/cns_candidates.txt"

echo "Working Folder: ${WRK_DIR}"
echo "Candidate Options: ${CAN_OPTIONS}"
echo "Consensus Options: ${CNS_OPTIONS}"
echo "Read List: ${READS_LIST}"

### pack database
if [ -f ${PACK_JOB_FINISHED} ]; then
    echo "Job [PackData] Has Been Finished, Skip It."
else
    pwd
    PACK_CMD="${PACK_WORKER} ${PACKED_DATA_DIR} ${READS_LIST}"
    OntCns2RunCmd.sh ${PACK_CMD}
    if [ $? -ne 0 ]; 
    then
        exit 1;
    fi
    touch ${PACK_JOB_FINISHED}
fi

### find candidates
if [ -f ${CAN_JOB_FINISHED} ]; then
    echo "Job [CanDetect] Has Been Finished, Skip It."
else
    CAN_CMD="${CAN_WORKER} ${CAN_OPTIONS} ${PACKED_DATA_DIR} ${CANDIDATES}"
    echo "CAN_CMD: ${CAN_CMD}"
    OntCns2RunCmd.sh ${CAN_CMD}
    if [ $? -ne 0 ];
    then 
        exit 1;
    fi
    touch ${CAN_JOB_FINISHED}
fi

### consensus
CNS_CMD="${CNS_WORKER} ${CNS_OPTIONS} ${PACKED_DATA_DIR} ${CANDIDATES} ${CNS_RESULTS} ${UNCNS_RESULTS}"
OntCns2RunCmd.sh ${CNS_CMD}
if [ $? -ne 0 ];
then
    exit 1;
fi
touch ${CNS_JOB_FINISHED}
