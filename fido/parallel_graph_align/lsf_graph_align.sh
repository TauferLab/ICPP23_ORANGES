#!/usr/bin/env bash 


n_runs=$1
n_procs=$2
input_graph1=$3
input_graph2=$4
results_path=$5

sims=1

function join_by { local d=$1; shift; local f=$1; shift; printf %s "$f" "${@/#/$d}"; }

# Define Paths
source ./fido_paths.config
mkdir -p ${results_path}
root_path=${fido_project_root}/fido
graph_align_job_script=${root_path}/parallel_graph_align/graph_align_par.sh
load_assignment="static"

# Configure Outputs
output_file=${results_path}/debugging/lsf_output.txt
error_file=${results_path}/debugging/lsf_error.txt

# Clean and Compile for New Run  
cd ${root_path}/parallel_graph_align
make clean
make all
make deg
echo "Made executable"


# Submit job to calculate graph alignment
echo ${results_path}
unq_job=$(date +%s)
bsub -o ${output_file} -e ${error_file} -n ${n_procs} -m "tellico-compute0" -J "fidojobs_${unq_job}[1-${sims}]" ${graph_align_job_script} ${n_runs} ${n_procs} ${input_graph1} ${input_graph2} ${load_assignment} ${results_path} -f input.\$LSB_JOBINDEX
bsub -n $((32-${n_procs}*${sims})) -m "tellico-compute0" -J "fidobuff" sleep 400000 
bwait -w "done(fidojobs_${unq_job})"
bkill -J "fidobuff"
echo "" >> ${root_path}/parallel_graph_align/time_results.txt

