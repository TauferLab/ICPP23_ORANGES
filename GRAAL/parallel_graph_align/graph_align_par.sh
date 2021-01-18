#!/usr/bin/env bash                                                                                                                                                                

#BSUB -o build_graph-%j.out                                                                                                                                                                                   
#BSUB -e build_graph-%j.err                                                                                                                                                                           

#n_procs=$1
#input_graph=$2
#graph_alignment_bin=$2
n_runs=$1
n_procs=$2
runs1=$3
runs2=$4
slice=$5
results_path=$6

root_path=$HOME/Src_GraphAlignment/GRAAL

# Define Paths for Input Files
#input_graph=${root_path}/test_graphs/graphn1.txt
#graph_folder=${root_path}/test_graphs
graph_folder=${root_path}/strong_scaling_txt
orbit_file=${root_path}/orbit_list.txt

# Define Standard IO
#std_out=${root_path}/debugging/graph_align_output.txt
#std_err=${root_path}/debugging/graph_align_error.txt

# Clean and Compile for New Run
#cd ${root_path}/parallel_graph_align
#make clean
#make all
#echo "Made executable"

#comm_procs=(10 16 24)
#comm_iters=(1 2 4 8 16)
#comm_msg_sizes=(512)

#comm_procs=(10 16 24)
#comm_iters=(4 2 1)
#comm_sizes=(92 92 94)

# Configure Inputs                     
#run_scales=(4 8 16 32 64)
#run_scales=(8)
#input_graphs_1=(slice_2.txt slice_3.txt slice_4.txt slice_5.txt slice_6.txt)
#input_graphs_2=(slice_2.txt slice_3.txt slice_4.txt slice_5.txt slice_6.txt)
input_graphs_1=(slice_2.txt slice_3.txt)
input_graphs_2=(slice_4.txt slice_5.txt)
#input_graphs_1=(sgraph1.txt sgraph2.txt)
#input_graphs_2=(sgraph1.txt sgraph2.txt)
#input_graphs=(graph_slice3.txt)
run_idx_low=1
run_idx_high=${n_runs}

echo "Set up Configuration Options" 
#for input_graph1 in ${input_graphs_1[@]};
#for input_graph1 in ${graph_folder}/*;
#do
    #for input_graph2 in ${input_graphs_2[@]};
    #for input_graph2 in ${graph_folder}/*;
    #do
#for procs in ${comm_procs[@]};
#do
    #for iters in ${comm_iters[@]};
    #do
	#${graph_folder}/slice2_message_race_niters_${iters}_nprocs_${procs}_msg_size_512_run_001_*.txt ${graph_folder}/slice2_message_race_niters_${iters}_nprocs_${procs}_msg_size_512_run_002_*.txt
	#echo "For Input Graph ${input_graph}"
#	for n_procs in ${run_scales[@]};
#	do
	    #echo "For ${n_procs} Processes"
        #for run_idx in `seq -f "%03g" ${run_idx_low} ${run_idx_high}`;
        #do
	    #echo "======================================================================"
	    #echo "Starting run ${run_idx} of ${input_graph1} with ${input_graph2} on ${n_procs} processes"
	    #mpirun -np ${n_procs} > ${std_out} 2> ${std_err} ./graph_alignment ${graph_folder}/slice2_message_race_niters_${iters}_nprocs_${procs}_msg_size_512_run_001_*.txt ${graph_folder}/slice2_message_race_niters_${iters}_nprocs_${procs}_msg_size_512_run_002_*.txt ${orbit_file}
	    #mpirun -np ${n_procs} > ${std_out} 2> ${std_err} ./graph_alignment ${graph_folder}/${input_graph1} ${graph_folder}/${input_graph2} ${orbit_file}
	    #echo "======================================================================"
	#done
    #done
    #    done
#done


#for i in ${!comm_procs[@]}; do
for run_idx in `seq -f "%03g" ${run_idx_low} ${run_idx_high}`; do

    run_path=${results_path}/procs_${n_procs}/slice_${slice}/comm_runs_${runs1}_w_${runs2}/run_${run_idx}/
    #mkdir ${results_path}/slice_${slice}/comm_runs_${runs1}_w_${runs2}/run_${run_idx}/
    mkdir -p ${run_path}
    cd ${run_path}

    echo "======================================================================"
    echo "Starting run ${run_idx} of ${input_graph1} with ${input_graph2} on ${n_procs} processes"
    mpirun -np ${n_procs} > ${run_path}/graph_align_out.txt 2> ${run_path}/graph_align_err.txt ${root_path}/parallel_graph_align/graph_alignment ${graph_folder}/slice_${slice}_message_race_niters_10_nprocs_16_msg_size_512_run_00${runs1}_size_332.txt ${graph_folder}/slice_${slice}_message_race_niters_10_nprocs_16_msg_size_512_run_00${runs2}_size_332.txt ${orbit_file} ${root_path}/parallel_graph_align/time_results.txt
    echo "======================================================================"

done
#done
