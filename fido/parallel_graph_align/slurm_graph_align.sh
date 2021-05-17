#!/usr/bin/env bash 


n_runs=$1
n_procs=$2
input_graph1=$3
input_graph2=$4
results_path=$5
sims=1
slurm_queue="skx-normal"
slurm_time_limit=$6

function join_by { local d=$1; shift; local f=$1; shift; printf %s "$f" "${@/#/$d}"; }

#echo "Entered file"
# Define Paths
source ./fido_paths.config
mkdir -p ${results_path}
root_path=${fido_project_root}/fido
graph_align_job_script=${root_path}/parallel_graph_align/graph_align_par.sh
load_assignment="static"
partition="metis"
#partition="divide"

# Configure Outputs
#echo "Configuring Outputs"
output_file=${results_path}/slurm_output.txt
error_file=${results_path}/slurm_error.txt

# Clean and Compile for New Run  
#echo "Changing directory"
cd ${root_path}/parallel_graph_align
#make clean
#make all
#make deg
#make stat
#echo "Made executable"


# Submit job to calculate graph alignment
echo ${results_path}
unq_job=$(date +%s)
sbatch --output ${output_file} --error ${error_file} -N 1 -J fido_job -p ${slurm_queue} -n ${n_procs} -t ${slurm_time_limit} ${graph_align_job_script} ${n_runs} ${n_procs} ${input_graph1} ${input_graph2} ${load_assignment} ${partition} ${results_path}
#sbatch -N 1 --wait -J fido_job -p ${slurm_queue} -n ${n_procs} -t ${slurm_time_limit} ${graph_align_job_script} ${n_runs} ${n_procs} ${input_graph1} ${input_graph2} ${load_assignment} ${partition} ${results_path}
echo "" >> ${root_path}/parallel_graph_align/time_results.txt

