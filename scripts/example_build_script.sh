#!/usr/bin/env bash

cmake \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DKokkos_ROOT=/home/ntan1/kokkos/build/install \
  -DKokkosKernels_DIR=/home/ntan1/kokkos-kernels/build/install/lib64/cmake/KokkosKernels \
  -Ddeduplicator_DIR=/home/ntan1/Src_Deduplication_Module/build/install/lib/cmake/deduplicator \
  ..

#  -DCMAKE_CXX_COMPILER=/home/ntan1/kokkos/bin/nvcc_wrapper \
#  -Dresilience_DIR=/home/ntan1/kokkos-resilience/build/install/lib64/cmake/resilience/ \
