#!/usr/bin/env bash                                                                                                                                                                

#BSUB -o build_graph-%j.out                                                                                                                                                                                   
#BSUB -e build_graph-%j.err                                                                                                                                                                           

n_procs=$1
graph_alignment_bin=$2

root_path=$HOME/Src_GraphAlignment/GRAAL

input_graph=${root_path}/graphn1.txt
orbit_file=${root_path}/orbit_list.txt

std_out=${root_path}/debugging/graph_align_output.txt
std_err=${root_path}/debugging/graph_align_error.txt

mpirun -np ${n_procs} ${graph_alignment_bin} ${input_graph} ${orbit_file} > ${std_out} 2> ${std_err}
