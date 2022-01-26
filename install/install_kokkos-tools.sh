#!/usr/bin/env bash


build_dir=build
#memory_usage_dir=$1

#rm -rf ${build_dir}
#mkdir -p ${build_dir}
#cd ${build_dir}

#outer_dir=$(pwd)
#cd ${memory_usage_dir}

#cmake .. \
#  -DCMAKE_INSTALL_PREFIX=$(pwd) \
#  -DKokkos_ROOT=${kokkos_build_dir}/${sys_lib}/cmake/Kokkos/
#make
#make install
chmod +x build-all.sh
./build-all.sh

#cd ${outer_dir}

