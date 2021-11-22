#!/usr/bin/env bash


build_dir=build/
project_root=$(pwd)
#kokkos_path=/work2/08092/tg873913/stampede2/spack/opt/spack/linux-centos7-skylake_avx512/intel-19.1.1.217/kokkos-3.3.01-4ct54q7mxt37js5622iokn63yrviwyjz/
#kokkos_kernels_path=/work2/08092/tg873913/stampede2/spack/opt/spack/linux-centos7-skylake_avx512/intel-19.1.1.217/kokkos-kernels-3.2.00-ooqvbeobrgycgsmt453f2v3lhqy7zuia/


## Setup progress delimiter
n_columns=$(stty size | awk '{print $2}')
progress_delimiter=""
for i in `seq 1 ${n_columns}`;
do
    progress_delimiter+="-"
done


## Update Git Submodules
echo
echo ${progress_delimiter}
echo "Fetching Submodules..."
echo ${progress_delimiter}
echo
rm -rf ./submodules/kokkos*
git submodule update --init --recursive
echo
echo ${progress_delimiter}
echo "Done Fetching Submodules."
echo ${progress_delimiter}
echo


## Link ESSENS
echo 
echo ${progress_delimiter}
echo "Linking ESSENS..."
echo ${progress_delimiter}
echo
cp fido/config/CMakeLists.txt submodules/ESSENS/
cp fido/config/EssensConfig.cmake.in submodules/ESSENS/
echo
echo ${progress_delimiter}
echo "Done Linking ESSENS"
echo ${progress_delimiter}
echo


## Build Kokkos
echo
echo ${progress_delimiter}
echo "Building Kokkos..."
echo ${progress_delimiter}
echo
kokkos_build_dir=${project_root}/submodules/kokkos/${build_dir}
cd submodules/kokkos/
. ${project_root}/install/install_kokkos.sh
cd ../../
echo
echo ${progress_delimiter}
echo "Done Installing Kokkos"
echo ${progress_delimiter}
echo


## Build Kokkos-Kernels
echo
echo ${progress_delimiter}
echo "Building Kokkos-Kernels..."
echo ${progress_delimiter}
echo
kokkos_kernels_build_dir=${project_root}/submodules/kokkos-kernels/${build_dir}
cd submodules/kokkos-kernels/
. ${project_root}/install/install_kokkos-kernels.sh ${kokkos_build_dir}
cd ../../
echo
echo ${progress_delimiter}
echo "Done Installing Kokkos-Kernels"
echo ${progress_delimiter}
echo


## Build fido
echo
echo ${progress_delimiter}
echo "Installing Fido..."
echo ${progress_delimiter}
echo
. ./install/install_fido.sh
echo
echo ${progress_delimiter}
echo "Done Installing Fido"
echo ${progress_delimiter}
echo


## Setup project directory
sed -i "s#fido_project_root= #fido_project_root=${project_root}#" ./fido/config/fido_paths.config



