#!/usr/bin/env bash


n_procs=$1
input_graph1=$2
input_graph2=$3
load_assignment=$4
paths_dir=$5
results_path=$6

#echo ${results_path}
#echo $#
#echo $@

source ${paths_dir}/fido_paths.config
#root_path=${fido_project_root}/fido
#g_align=${root_path}/parallel_graph_align/graph_alignment
#g_align=${fido_project_root}/build/Fido
#deg_count=${root_path}/parallel_graph_align/node_deg

# Define Paths for Input Files
#orbit_file=${root_path}/orbit_list.txt
#time_keeping=${root_path}/parallel_graph_align/time_results.txt

#val_mpi_suppr=/home/pnbell/spack/opt/spack/linux-rhel7-power9le/gcc-9.2.0/openmpi-3.1.6-3kn5q2dl3x5lk4jhkni75o4zx5dijrfq/share/openmpi/openmpi-valgrind.supp
#val_mpi_suppr2=/home/pnbell/spack/opt/spack/linux-rhel7-power9le/gcc-9.2.0/openmpi-3.1.6-klg2sl3myyo6etpxcj4os46qhd2bahfh/share/openmpi/openmpi-valgrind.supp

#echo ${input_graph1}
#echo ${input_graph2}

#run_idx_low=0
#run_idx_high=$((n_runs-1))

export PSM2_MEMORY=large

#for run_idx in `seq -f "%03g" ${run_idx_low} ${run_idx_high}`; do
    
# Prepare output path
#echo ${results_path}
#run_path=${results_path}/sims_$LSB_JOBINDEX/procs_${n_procs}/run_${run_idx}/
run_path=${results_path}/procs_${n_procs}/
mkdir -p ${run_path}
cd ${run_path}
#mkdir -p ${run_path}/runtime_data

debug_path=${run_path}/debug/
mkdir -p ${debug_path}

runtime_path=${run_path}/runtime_data/
mkdir -p ${runtime_path}

echo "======================================================================"
echo "Starting Fido alignment of ${input_graph1} with ${input_graph2} on ${n_procs} processes"
    	#mpirun -np ${n_procs} > ${run_path}/graph_align_out.txt 2> ${run_path}/graph_align_err.txt valgrind --leak-check=full --error-limit=no --log-file="valgrind_out.txt" --suppressions=${val_mpi_suppr2} --suppressions=${val_mpi_suppr} --suppressions=/home/pnbell/Src_GraphAlignment/GRAAL/parallel_graph_align/mpi_supp_samp.supp --gen-suppressions=all ${g_align} ${input_graph1} ${input_graph2} ${orbit_file} ${time_keeping}
	#mpirun -np ${n_procs} > ${run_path}/graph_align_out.txt 2> ${run_path}/graph_align_err.txt valgrind --leak-check=full --error-limit=no --log-file="valgrind_out.txt" --suppressions=${val_mpi_suppr2} --suppressions=${val_mpi_suppr} --suppressions=/home/pnbell/Src_GraphAlignment/GRAAL/parallel_graph_align/mpi_supp_samp.supp ${g_align} ${input_graph1} ${input_graph2} ${orbit_file} ${time_keeping}
mpirun -np ${n_procs} > ${debug_path}/graph_align_out.txt 2> ${debug_path}/graph_align_err.txt ${g_align} ${input_graph1} ${input_graph2} ${orbit_file} --kokkos-threads=1
	#mpirun -np ${n_procs} > ${run_path}/deg_count_out.txt 2> ${run_path}/deg_count_err.txt ${deg_count} ${input_graph1} ${input_graph2}
	#mpirun -np ${n_procs} > graph_align_out.txt 2> graph_align_err.txt ${g_align} ${input_graph1} ${input_graph2} ${orbit_file} ${time_keeping}
echo "======================================================================"

#done

# Construct deg distribution across processes
if [ "${load_assignment}" == "static" ]; then
    cd ${run_path}/../
    mpirun -np ${n_procs} > deg_count_out.txt 2> deg_count_err.txt ${deg_count} ${input_graph1} ${input_graph2}
fi
