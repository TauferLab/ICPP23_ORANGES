#!/usr/bin/env bash


build_dir=build/
project_root=$(pwd)

kokkos_path=/work2/08092/tg873913/stampede2/spack/opt/spack/linux-centos7-skylake_avx512/intel-19.1.1.217/kokkos-3.3.01-4ct54q7mxt37js5622iokn63yrviwyjz/
kokkos_kernels_path=/work2/08092/tg873913/stampede2/spack/opt/spack/linux-centos7-skylake_avx512/intel-19.1.1.217/kokkos-kernels-3.2.00-ooqvbeobrgycgsmt453f2v3lhqy7zuia/


## Setup progress delimiter
n_columns=$(stty size | awk '{print $2}')
progress_delimiter=""
for i in `seq 1 ${n_columns}`;
do
    progress_delimiter+="-"
done


## Build fido
echo
echo ${progress_delimiter}
echo "Installing Fido"
echo ${progress_delimiter}
echo
rm -rf ${build_dir} && mkdir ${build_dir}
cd ${build_dir}
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
cd ..
echo
echo ${progress_delimiter}
echo "Done Installing Fido"
echo ${progress_delimiter}
echo


## Setup project directory
sed -i "s#fido_project_root= #fido_project_root=${project_root}#" ./fido/parallel_graph_align/fido_paths.config
