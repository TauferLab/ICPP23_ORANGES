#!/usr/bin/env bash


while [ -n "$1" ]; do
        case "$1" in
		-r) track_runtime_cmake=$2; shift; shift ;;
		#-t) num_threads_cmake=$2; shift; shift ;;
		-s) sim_mat_cutoff=$2; shift; shift ;;
		-l) library_type=$2; shift; shift ;;
		-c) num_chunks=$2; shift; shift ;;
	esac
done

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


## Process User Input
#echo 
#echo ${progress_delimiter}
#echo "Processing User Input"
#echo ${progress_delimiter}
#echo
#if [[ "${track_runtime}" == "True" ]]; then
#        track_runtime_cmake=1
#else
#        track_runtime_cmake=0
#fi
#echo
#echo ${progress_delimiter}
#echo "Done Processing User Input"
#echo ${progress_delimiter}
#echo


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
. ${project_root}/install/install_kokkos-kernels.sh ${kokkos_build_dir} ${library_type}
cd ../../
echo
echo ${progress_delimiter}
echo "Done Installing Kokkos-Kernels"
echo ${progress_delimiter}
echo


echo
echo ${progress_delimiter}
echo "Building Kokkos-Tools..."
echo ${progress_delimiter}
echo
kokkos_tools_dir=${project_root}/submodules/kokkos-tools/
#memory_usage_dir=${kokkos_tools_dir}
cd submodules/kokkos-tools/
. ${project_root}/install/install_kokkos-tools.sh
cd ../../
echo
echo ${progress_delimiter}
echo "Done Installing Kokkos-Tools"
echo ${progress_delimiter}
echo


## Build fido
echo
echo ${progress_delimiter}
echo "Installing Fido..."
echo ${progress_delimiter}
echo
. ./install/install_fido.sh ${track_runtime_cmake} ${sim_mat_cutoff} ${num_chunks} ${library_type}
echo
echo ${progress_delimiter}
echo "Done Installing Fido"
echo ${progress_delimiter}
echo


## Setup project directory
sed -i "s#fido_project_root= #fido_project_root=${project_root}#" ./fido/config/fido_paths.config
#sed -i "s#kokkos_tools_memory_usage= #kokkos_tools_memory_usage=${kokkos_tools_dir}#" ./fido/config/fido_paths.config


