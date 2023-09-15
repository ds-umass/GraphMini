#!/bin/bash
CURRENT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BUILD_TYPE=Release
# Build GraphMini
PROJECT_NAME="GraphMini"
mkdir -p ${CURRENT_DIR}/${PROJECT_NAME}/build
rm -rf ${CURRENT_DIR}/${PROJECT_NAME}/build/*
cmake -S ${CURRENT_DIR}/${PROJECT_NAME} -B ${CURRENT_DIR}/${PROJECT_NAME}/build -DCMAKE_BUILD_TYPE=${BUILD_TYPE} 
cmake --build ${CURRENT_DIR}/${PROJECT_NAME}/build -j;

# Build GraphPi
PROJECT_NAME="GraphPi"
mkdir -p ${CURRENT_DIR}/${PROJECT_NAME}/build
rm -rf ${CURRENT_DIR}/${PROJECT_NAME}/build/*
cmake -S ${CURRENT_DIR}/${PROJECT_NAME} -B ${CURRENT_DIR}/${PROJECT_NAME}/build -DCMAKE_BUILD_TYPE=${BUILD_TYPE} 
cmake --build ${CURRENT_DIR}/${PROJECT_NAME}/build -j;

# Build Dryadic
cd ${CURRENT_DIR}/Dryadic && make;
cd ${CURRENT_DIR}/Dryadic/graph_converter && make -j;
cd ${CURRENT_DIR}