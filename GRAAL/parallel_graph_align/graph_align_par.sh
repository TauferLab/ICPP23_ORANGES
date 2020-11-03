#!/usr/bin/env bash                                                                                                                                                                

#BSUB -o build_graph-%j.out                                                                                                                                                                                   
#BSUB -e build_graph-%j.err                                                                                                                                                                           

#n_procs=$1
#input_graph=$2
#graph_alignment_bin=$2
n_runs=$1

root_path=$HOME/Src_GraphAlignment/GRAAL

# Define Paths for Input Files
#input_graph=${root_path}/test_graphs/graphn1.txt
graph_folder=${root_path}/test_graphs
orbit_file=${root_path}/orbit_list.txt

# Define Standard IO
std_out=${root_path}/debugging/graph_align_output.txt
std_err=${root_path}/debugging/graph_align_error.txt

# Clean and Compile for New Run
cd ${root_path}/parallel_graph_align
make clean
make all
echo "Made executable"

# Configure Inputs                     
#run_scales=(4 8 16 32 64)
run_scales=(2)
input_graphs=(graphn1.txt)
#input_graphs=(graph_slice3.txt)
run_idx_low=1
run_idx_high=${n_runs}

echo "Set up Configuration Options" 
for input_graph in ${input_graphs[@]};
do
    #echo "For Input Graph ${input_graph}"
    for n_procs in ${run_scales[@]};
    do
	#echo "For ${n_procs} Processes"
        for run_idx in `seq -f "%03g" ${run_idx_low} ${run_idx_high}`;
        do
	    echo "======================================================================"
	    echo "Starting run ${run_idx} of ${input_graph} on ${n_procs} processes"
	    mpirun --oversubscribe -np ${n_procs} > ${std_out} 2> ${std_err} ./graph_alignment ${graph_folder}/${input_graph} ${orbit_file}
	    echo "======================================================================"
	done
    done
done
