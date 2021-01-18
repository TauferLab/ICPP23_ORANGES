#!/usr/bin/env bash 

#n_procs=$1
#input_graph=$1
n_runs=$1
results_path=$2

# Define Paths
root_path=$HOME/Src_GraphAlignment/GRAAL
graph_align_job_script=${root_path}/parallel_graph_align/graph_align_par.sh

# Configure Outputs
output_file=${root_path}/debugging/lsf_output.txt
error_file=${root_path}/debugging/lsf_error.txt

# Clean and Compile for New Run  
cd ${root_path}/parallel_graph_align
make clean
make all
echo "Made executable"

# Configure Inputs
#run_scales=(4 8 16 32 64)
#run_scales=(2 4)
#run_scales=(32)
#run_scales=(64)
#run_scales=(16 32 64)
#run_scales=(24 40 48 56)
#run_scales=(1 2 4 8)
#run_scales=(1 2 4 8 16 32 64)
run_scales=(1 2 4)

comm_procs=16
comm_iters=10
comm_size=332
comm_runs=3
comm_slices=2

#input_graphs=(graph_slice3.txt)
#run_idx_low=1
#run_idx_high=${n_runs}


#for run_idx in `seq -f "%03g" ${run_idx_low} ${run_idx_high}`; 
#do
#for input_graph in ${input_graph[@]};
#do
for n_procs in ${run_scales[@]}; do
#	for run_idx in `seq -f "%03g" ${run_idx_low} ${run_idx_high}`;
#	do
    for runs1 in $(seq 1 ${comm_runs}); do
	for runs2 in $(seq 1 ${comm_runs}); do
	    if [ ${runs1} -lt ${runs2} ]; then
		#echo ${runs1}, ${runs2}
		for slices in $(seq 1 ${comm_slices}); do
		    bsub -o ${output_file} -e ${error_file} -n ${n_procs} ${graph_align_job_script} ${n_runs} ${n_procs} ${runs1} ${runs2} ${slices} ${results_path}
		done
	    fi
	done
    done
done
