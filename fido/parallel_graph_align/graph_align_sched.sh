#!/usr/bin/env bash


n_procs=$1
n_nodes=$2
input_graph1=$3
input_graph2=$4
paths_dir=$5
results_path=$6
#sims=1
scheduler=$7
sched_queue=$8
sched_time_limit=$9
load_assignment=${10}

function join_by { local d=$1; shift; local f=$1; shift; printf %s "$f" "${@/#/$d}"; }

# Define Paths
source ${paths_dir}/fido_paths.config
mkdir -p ${results_path}
#root_path=${fido_project_root}/fido
#graph_align_job_script=${root_path}/parallel_graph_align/graph_align_par.sh
#load_assignment="dynamic"

# Configure Outputs
output_file=${results_path}/sched_output.txt
error_file=${results_path}/sched_error.txt

# Clean and Compile for New Run  
cd ${root_path}/parallel_graph_align
#make clean
#make all
#make deg
#echo "Made executable"


# Submit job to calculate graph alignment
#echo ${results_path}
unq_job=$(date +%s)
n_procs_per_node=$((n_procs/n_nodes))

if [ "${scheduler}" == "slurm" ]; then
	sbatch -o ${output_file} -e ${error_file} -N ${n_nodes} --ntasks-per-node=$((n_procs_per_node+1)) --wait -J fido_job -p ${sched_queue} -n ${n_procs} -t ${sched_time_limit} ${graph_align_job_script} ${n_procs} ${input_graph1} ${input_graph2} ${load_assignment} ${paths_dir} ${results_path}
fi
if [ "${scheduler}" == "lsf" ]; then
	bsub -o ${output_file} -e ${error_file} -n ${n_procs} -J "fidojobs_${unq_job}" -R "span[ptile=$((n_procs_per_node+1))]" -q ${sched_queue} -W ${sched_time_limit} ${graph_align_job_script} ${n_procs} ${input_graph1} ${input_graph2} ${load_assignment} ${paths_dir} ${results_path}
fi
if [ "${scheduler}" == "none" ]; then
	bash > ${output_file} 2> ${error_file} ${graph_align_job_script} ${n_procs} ${input_graph1} ${input_graph2} ${load_assignment} ${paths_dir} ${results_path}
fi

#echo "" >> ${root_path}/parallel_graph_align/time_results.txt
