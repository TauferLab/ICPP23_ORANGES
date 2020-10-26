#!/usr/bin/env bash 

n_procs=$1
input_graph=$2

root_path=$HOME/Src_GraphAlignment/GRAAL

graph_align_job_script=${root_path}/parallel_graph_align/graph_align_par_lsf.sh

output_file=${root_path}/debugging/lsf_output.txt
error_file=${root_path}/debugging/lsf_error.txt

bsub > ${output_file} 2> ${error_file} ${graph_align_job_script} ${n_procs} ${input_graph}
