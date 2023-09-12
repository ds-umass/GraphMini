#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

bash ${SCRIPT_DIR}/download.sh
bash ${SCRIPT_DIR}/graphpi_prep.sh
bash ${SCRIPT_DIR}/minigraph_prep.sh
bash ${SCRIPT_DIR}/dryadic_prep.sh
