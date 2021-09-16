#!/usr/bin/env bash


kokkos_path=$1
kokkos_kernels_path=$2

export CC=mpicc
export CXX=mpicxx
cmake \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_PREFIX_PATH=${project_root} \
  -DKokkos_DIR=${kokkos_path}/lib64/cmake/Kokkos \
  -DKokkosKernels_DIR=${kokkos_kernels_path}/lib64/cmake/KokkosKernels \
  ..
make
