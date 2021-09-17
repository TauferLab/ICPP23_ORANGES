#!/usr/bin/env bash

build_dir=build

rm -rf ${build_dir}
mkdir -p ${build_dir}
cd ${build_dir}

export CC=mpicc
export CXX=mpicxx
cmake \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_INSTALL_PREFIX=$(pwd) \
  -DKokkos_ENABLE_OPENMP=ON \
  -DKokkos_ENABLE_SERIAL=ON \
  .. 
make
make install

cd ../
