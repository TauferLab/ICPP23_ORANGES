#!/usr/bin/env bash 

#n_procs=$1
#input_graph=$1
n_runs=$1

# Define Paths
root_path=$HOME/Src_GraphAlignment/GRAAL
graph_align_job_script=${root_path}/parallel_graph_align/graph_align_par.sh

# Configure Outputs
output_file=${root_path}/debugging/lsf_output.txt
error_file=${root_path}/debugging/lsf_error.txt

# Configure Inputs
run_scales=(4 8 16 32 64)
input_graphs=(graph_slice3.txt)
run_idx_low=1
run_idx_high=${n_runs}

#for run_idx in `seq -f "%03g" ${run_idx_low} ${run_idx_high}`; 
#do
for input_graph in ${input_graph[@]};
do
    for n_procs in ${run_scales[@]};
    do
	for run_idx in `seq -f "%03g" ${run_idx_low} ${run_idx_high}`;
	do
	    bsub -o ${output_file} -e ${error_file} -n ${n_procs} ${graph_align_job_script} ${n_procs} ${input_graph}
	done
    done
done
