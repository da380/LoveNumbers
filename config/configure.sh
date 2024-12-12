#!/bin/bash

BUILD_EXAMPLES="ON"

MPI_PATH=""
MFEM_DIR=""

cmake -S . -B build \
           -DMFEM_DIR=$MFEM_DIR \
           -DMPI_C_COMPILER="${MPI_PATH}mpicc" \
           -DMPI_CXX_COMPILER="${MPI_PATH}mpic++" \
           -DBUILD_EXAMPLES=$BUILD_EXAMPLES

