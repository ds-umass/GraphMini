#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
DATA_INDIR=${SCRIPT_DIR}/../Datasets/GraphMini
EXPERIMENT_OUTDIR=${SCRIPT_DIR}/Base
BIN_PATH=${SCRIPT_DIR}/../GraphMini/build/bin/run
source ${SCRIPT_DIR}/queries.sh

adj_type=1
prun_type=0
par_type=0
QUERY_TYPE="EdgeInduced"

for graph_name in "${GraphNames[@]}"
    do 
        INDEX=1
        for query in "${Queries[@]}"
            do
                start_time=`date +%Y-%m-%d-%T`
                query_name="P${INDEX}"
                query_size=${QuerySizes[${INDEX} - 1]}
                path_to_graph=${DATA_INDIR}/${graph_name}/
                output_dir=${EXPERIMENT_OUTDIR}/${graph_name}/${query_name}/${QUERY_TYPE}
                output_path=${output_dir}/run.log
                mkdir -p ${output_dir}

                start_time=`date +%Y-%m-%d-%T`
                start_time_ms=$(date +%s.%N)
                echo "--------- START EXPIERMENT ${start_time} ----------" >> ${output_path}
                
                timeout ${TIMEOUT} ${BIN_PATH} ${graph_name} ${path_to_graph} ${query_name} ${query} ${adj_type} ${prun_type} ${par_type} >> ${output_path}

                end_time_ms=$(date +%s.%N)
                end_time=`date +%Y-%m-%d-%T`                
                duration=$(echo "($end_time_ms - $start_time_ms)" | bc)
                echo "--------- END EXPIERMENT ${end_time} (${duration} seconds)----------" >> ${output_path}
                INDEX=$((INDEX+1))
            done
    done