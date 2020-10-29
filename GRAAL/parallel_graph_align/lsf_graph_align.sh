#!/usr/bin/env bash 

#n_procs=$1
input_graph=$1

root_path=$HOME/Src_GraphAlignment/GRAAL

graph_align_job_script=${root_path}/parallel_graph_align/graph_align_par.sh

output_file=${root_path}/debugging/lsf_output.txt
error_file=${root_path}/debugging/lsf_error.txt

run_scales=(32)

for n_procs in ${run_scales[@]};
do
    bsub -o ${output_file} -e ${error_file} -n ${n_procs} ${graph_align_job_script} ${n_procs} ${input_graph}
done
