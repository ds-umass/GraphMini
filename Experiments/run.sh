#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

bash ${SCRIPT_DIR}/run_graphmini_edge_iep.sh
bash ${SCRIPT_DIR}/run_graphmini_edge.sh
bash ${SCRIPT_DIR}/run_graphmini_vertex.sh

bash ${SCRIPT_DIR}/run_graphpi_edge_iep.sh
bash ${SCRIPT_DIR}/run_graphpi_edge.sh

bash ${SCRIPT_DIR}/run_dryadic_edge.sh
bash ${SCRIPT_DIR}/run_dryadic_vertex.sh

bash ${SCRIPT_DIR}/run_base_edge_iep.sh
bash ${SCRIPT_DIR}/run_base_edge.sh
bash ${SCRIPT_DIR}/run_base_vertex.sh