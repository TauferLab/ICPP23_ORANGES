#!/usr/bin/env bash                                                                                                                                                                

#BSUB -o build_graph-%j.out                                                                                                                                                                                   
#BSUB -e build_graph-%j.err                                                                                                                                                                           

n_procs=$1
input_graph=$2
#graph_alignment_bin=$2

root_path=$HOME/Src_GraphAlignment/GRAAL

# Define Paths for Input Files
#input_graph=${root_path}/test_graphs/graphn1.txt
graph_folder=${root_path}/test_graphs
orbit_file=${root_path}/orbit_list.txt

# Define Standard IO
std_out=${root_path}/debugging/graph_align_output.txt
std_err=${root_path}/debugging/graph_align_error.txt

# Clean and Compile for New Run
make clean
make all

mpirun -np ${n_procs} ./graph_alignment ${graph_folder}/${input_graph} ${orbit_file} > ${std_out} 2> ${std_err}
