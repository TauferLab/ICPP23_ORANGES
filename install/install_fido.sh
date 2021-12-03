#!/usr/bin/env bash

build_dir=build
sys_lib=lib
project_root=$(pwd)

rm -rf ${build_dir}
mkdir -p ${build_dir}
cd ${build_dir}

#kokkos_path=$1
#kokkos_kernels_path=$2

export CC=mpicc
export CXX=mpicxx
cmake \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DTRACK_RUNTIME:STRING=0 \
  -DTHREAD_COUNT:STRING=1 \
  -DCMAKE_PREFIX_PATH=${project_root} \
  -DKokkos_DIR=${project_root}/submodules/kokkos/${build_dir}/${sys_lib}/cmake/Kokkos \
  -DKokkosKernels_DIR=${project_root}/submodules/kokkos-kernels/${build_dir}/${sys_lib}/cmake/KokkosKernels \
  ..
make

cd ..
