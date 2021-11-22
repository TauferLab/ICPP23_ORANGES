#!/usr/bin/env bash


build_dir=build
kokkos_build_dir=$1
sys_lib=lib

rm -rf ${build_dir}
mkdir -p ${build_dir}
cd ${build_dir}

cmake .. \
  -DCMAKE_INSTALL_PREFIX=$(pwd) \
  -DKokkos_ROOT=${kokkos_build_dir}/${sys_lib}/cmake/Kokkos/
make
make install

cd ../

