#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BIN_PATH=${SCRIPT_DIR}/bin/query_l_best
GraphPath=$1
${BIN_PATH} ${GraphPath}