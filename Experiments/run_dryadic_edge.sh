#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source ${SCRIPT_DIR}/queries.sh

DATA_INDIR=${SCRIPT_DIR}/../Datasets/Dryadic
EXPERIMENT_OUTDIR=${SCRIPT_DIR}/Dryadic
DRYADIC_HOME=${SCRIPT_DIR}/../Dryadic
COMPILER_PATH=${DRYADIC_HOME}/compile.sh
RUNNT_PATH=${DRYADIC_HOME}/run.sh

QUERY_TYPE_FLAG=0
QUERY_TYPE="EdgeInduced"
TIMEOUT=24h
for graph_name in "${GraphNames[@]}"
    do 
        INDEX=1
        for query in "${Queries[@]}"
            do
                query_name="P${INDEX}"
                query_size=${QuerySizes[${INDEX} - 1]}
                path_to_graph=${DATA_INDIR}/${graph_name}/snap.txt
                output_dir=${EXPERIMENT_OUTDIR}/${graph_name}/${query_name}/${QUERY_TYPE}
                output_path=${output_dir}/run.log
                mkdir -p ${output_dir}

                start_time=`date +%Y-%m-%d-%T`
                start_time_ms=$(date +%s.%N)

                echo "--------- START EXPIERMENT ${start_time} ----------" >> ${output_path}
                
                timeout ${TIMEOUT} bash ${COMPILER_PATH} ${query_size} ${query} ${QUERY_TYPE_FLAG} >> ${output_path}
                timeout ${TIMEOUT} bash ${RUNNT_PATH} ${path_to_graph} >> ${output_path}
                
                end_time_ms=$(date +%s.%N)
                end_time=`date +%Y-%m-%d-%T`                
                duration=$(echo "($end_time_ms - $start_time_ms)" | bc)
                
                echo "--------- END EXPIERMENT ${end_time} (${duration} seconds)----------" >> ${output_path}
                INDEX=$((INDEX+1))
            done
    done