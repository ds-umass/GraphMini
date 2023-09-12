#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
QuerySize=$1
QueryStr=$2
QueryType=$3 # 0 : EdgeInduced; 1: VertexInduced

# generating code
QueryFile=${SCRIPT_DIR}/compiler/query
python3 ${SCRIPT_DIR}/get_query.py $QuerySize $QueryStr > ${QueryFile}
COMPILER_DIR=${SCRIPT_DIR}/compiler/

if [ "$QueryType" -eq 0 ]; then
    # echo "Query is edge induced"
    echo "#define EDGE_INDUCED" > ${SCRIPT_DIR}/query_t.h
else
    # echo "Query is vertex induced"
    echo "#undef EDGE_INDUCED" > ${SCRIPT_DIR}/query_t.h
fi

# compile to binary
start_time=$(date +%s.%N)
cd $COMPILER_DIR && ./compile --pattern-list-l ./query 1>${SCRIPT_DIR}/log/compile.log 2>${SCRIPT_DIR}/log/compile_err.log
cd ${SCRIPT_DIR} && make 1> ${SCRIPT_DIR}/log/make.log 2>${SCRIPT_DIR}/log/make_err.log
end_time=$(date +%s.%N)
duration=$(echo "($end_time - $start_time)" | bc)
echo "CODE_GENERATION_TIME(s)=$duration"