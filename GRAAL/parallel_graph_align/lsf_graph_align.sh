#!/usr/bin/env bash 

#n_procs=$1
#input_graph=$1
n_runs=$1
results_path=$2

function join_by { local d=$1; shift; local f=$1; shift; printf %s "$f" "${@/#/$d}"; }

# Define Paths
root_path=$HOME/Src_GraphAlignment/GRAAL
graph_align_job_script=${root_path}/parallel_graph_align/graph_align_par.sh

# Configure Outputs
output_file=${root_path}/debugging/lsf_output.txt
error_file=${root_path}/debugging/lsf_error.txt

# Clean and Compile for New Run  
#rm ${output_file}
#rm ${error_file}
#rm -r ${results_path}/*

cd ${root_path}/parallel_graph_align
make clean
make all
make deg
echo "Made executable"

# Configure Inputs
#run_scales=(4 8 16 32 64)
#run_scales=(2 4)
#run_scales=(32)
#run_scales=(64)
#run_scales=(16 32 64)
#run_scales=(24 40 48 56)
#run_scales=(1 8 16)
#run_scales=(1 2 4 8 16 32)
run_scales=(16 32)
#run_scales=(1)

sim_jobs=(1)
#sim_jobs=(16)
#sim_jobs=(1)
#sim_jobs=(32)

comm_type=amg2013
comm_procs=16
comm_iters=2
comm_size=1952
msg_size=1024
comm_runs=2
comm_slices=1

#input_graphs=(graph_slice3.txt)
#run_idx_low=1
#run_idx_high=${n_runs}

#job_sub=$( bsub -o ${output_file} -e ${error_file} -n 1 -m "tellico-compute1" ${graph_align_job_script} ${n_runs} ${n_procs} ${runs1} ${runs2} ${slices} ${iters} ${results_path} )
#job_sub_id=$( echo ${job_sub} | sed 's/[^0-9]*//g' )
#job_sub_arr+=("done(${job_sub})")
#job_join_str=$( join_by "&&" ${job_sub_arr[@]} )

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
		    for sims in ${sim_jobs[@]}; do
			#for iter in $(seq 0 $((sims-1))); do
			#for iters in `seq -f "%03g" 1 32`; do
			echo ${results_path}
			unq_job=$(date +%s)
		    	bsub -o ${output_file} -e ${error_file} -n ${n_procs} -m "tellico-compute0" -J "fidojobs_${unq_job}[1-${sims}]" ${graph_align_job_script} ${n_runs} ${n_procs} ${runs1} ${runs2} ${slices} ${comm_procs} ${comm_iters} ${comm_size} ${comm_type} ${msg_size} ${sims} ${results_path} -f input.\$LSB_JOBINDEX
			bsub -n $((32-${n_procs}*${sims})) -m "tellico-compute0" -J "fidobuff" sleep 400000 
			bwait -w "done(fidojobs_${unq_job})"
			bkill -J "fidobuff"
		    	    #job_sub_id=$( echo ${job_sub} | sed 's/[^0-9]*//g' )
			    #job_sub_arr+=("done(${job_sub})")
			    #job_join_str=$( join_by "&&" ${job_sub_arr[@]} )
			#done
			#done
			#wait
			echo "" >> ${root_path}/parallel_graph_align/time_results.txt
		    done
		    #echo "" >> ${root_path}/parallel_graph_align/time_results.txt 
		done
	    fi
	done
    done
done
