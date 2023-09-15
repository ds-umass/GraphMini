#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
DATASET_INDIR=${SCRIPT_DIR}/TXT
DATASET_OUTDIR=${SCRIPT_DIR}/GraphMini
BIN_PATH=${SCRIPT_DIR}/../GraphMini/build/bin/prep

echo "Preprocess graph files for GraphMini"
echo "The preprocessed files are written into" ${DATASET_OUTDIR}

graph_name="wiki"

for graph_name in "wiki" "youtube" "patents" "orkut" "lj" "friendster"
    do
        if test -f "${DATASET_OUTDIR}/${graph_name}/snap.txt"; then
            echo "${graph_name} exists"
        else
            mkdir -p ${DATASET_OUTDIR}/${graph_name}
            ln -s ${DATASET_INDIR}/${graph_name}/snap.txt ${DATASET_OUTDIR}/${graph_name}/snap.txt
            ${BIN_PATH} ${DATASET_OUTDIR}/${graph_name}/
        fi
    done