#!/usr/bin/env bash 

n_runs=$1
n_procs=$2
input_graph1=$3
input_graph2=$4
results_path=$5

function join_by { local d=$1; shift; local f=$1; shift; printf %s "$f" "${@/#/$d}"; }

# Define Paths
source ./fido_paths.config
mkdir -p ${results_path}
root_path=${fido_project_root}/fido
graph_align_job_script=${root_path}/parallel_graph_align/graph_align_par.sh
load_assignment="static"

# Clean and Compile for New Run
cd ${root_path}/parallel_graph_align
make clean
make all
make deg
echo "Made executable"


# Submit job to calculate graph alignment
echo ${results_path}
unq_job=$(date +%s)
bash ${graph_align_job_script} ${n_runs} ${n_procs} ${input_graph1} ${input_graph2} ${load_assignment} ${results_path}
echo "" >> ${root_path}/parallel_graph_align/time_results.txt

