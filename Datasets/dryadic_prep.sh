#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
DATASET_INDIR=${SCRIPT_DIR}/TXT
DATASET_OUTDIR=${SCRIPT_DIR}/Dryadic
BIN_DIR=${SCRIPT_DIR}/../Dryadic/graph_converter/

echo "Preprocess graph files for Dryadic"
echo "The preprocessed files are written into" ${DATASET_OUTDIR}

graph_name="wiki"

for graph_name in "wiki" "youtube" "patents" "orkut" "lj" "friendster"
    do
        if test -f "${DATASET_OUTDIR}/${graph_name}/snap.txt"; then
            echo "wiki exists"
        else
            mkdir -p ${DATASET_OUTDIR}/${graph_name}
            ln -s ${DATASET_INDIR}/${graph_name}/snap.txt ${DATASET_OUTDIR}/${graph_name}/snap.txt
            cd ${BIN_DIR} && bash ./label_graph.sh ${DATASET_OUTDIR}/${graph_name}/
        fi
    done