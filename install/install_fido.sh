#!/usr/bin/env bash

track_runtime=$1
num_threads=$2

build_dir=build
sys_lib=lib
project_root=$(pwd)

rm -rf ${build_dir}
mkdir -p ${build_dir}
cd ${build_dir}

#kokkos_path=$1
#kokkos_kernels_path=$2

track_runtime_cmake="${track_runtime_cmake:="1"}"
num_threads_cmake="${num_threads:="1"}"

export CC=mpicc
export CXX=mpicxx
cmake \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DTRACK_RUNTIME:STRING=${track_runtime_cmake} \
  -DTHREAD_COUNT:STRING=${num_threads_cmake} \
  -DCMAKE_PREFIX_PATH=${project_root} \
  -DKokkos_DIR=${project_root}/submodules/kokkos/${build_dir}/${sys_lib}/cmake/Kokkos \
  -DKokkosKernels_DIR=${project_root}/submodules/kokkos-kernels/${build_dir}/${sys_lib}/cmake/KokkosKernels \
  ..
make

cd ..
