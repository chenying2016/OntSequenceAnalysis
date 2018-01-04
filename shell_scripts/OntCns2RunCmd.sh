#!/bin/bash

i=0;
for O in $*;
do
    if [ $i -eq 0 ]; then
        CMD=$O
    else
        CMD="${CMD} $O"
    fi
    let i=i+1
done

echo "===$(date)=== [${CMD}] BEGINS"
${CMD}
if [ $? -ne 0 ]; then
    echo Failed At Running [${CMD}]
    exit 1;
fi
echo "===$(date)=== [${CMD}] FINISH"
