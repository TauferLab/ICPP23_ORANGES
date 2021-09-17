#!/usr/bin/env bash


build_dir=build
kokkos_build_dir=$1

rm -rf ${build_dir}
mkdir -p ${build_dir}
cd ${build_dir}

cmake .. \
  -DCMAKE_INSTALL_PREFIX=$(pwd) \
  -DKokkos_ROOT=${kokkos_build_dir}/lib64/cmake/Kokkos/
make
make install

cd ../

