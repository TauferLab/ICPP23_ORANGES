#!/usr/bin/env bash                                                                                                                                                                


n_runs=$1
n_procs=$2
input_graph1=$3
input_graph2=$4
load_assignment=$5
partition=$6
results_path=$7


#partition="divide"

source ./fido_paths.config
mkdir -p ${graphml_graph_files}
mkdir -p ${metis_graph_files}
mkdir -p ${part_dir}
graph_base_1=$(basename ${input_graph1} .txt)
graph_base_2=$(basename ${input_graph2} .txt)

run_idx_low=0
run_idx_high=$((n_runs-1))

# Manage Graph Paritioning
if [ ${partition} == "metis" ]; then
    
    # Make a graphml file for intermediary conversion
    ipython3 ${igraph_conversion_script} ${input_graph1} ${graphml_graph_files}
    ipython3 ${igraph_conversion_script} ${input_graph2} ${graphml_graph_files}
    metis_graphml1=${graphml_graph_files}/${graph_base_1}.graphml
    metis_graphml2=${graphml_graph_files}/${graph_base_2}.graphml

    # Make a .graph file for METIS
    ipython3 ${nk_conversion_script} ${metis_graphml1} ${metis_graph_files}
    #mv ${metis_graph_files}/event_graph.graph ${metis_graph_files}/event_graph_1.graph
    ipython3 ${nk_conversion_script} ${metis_graphml2} ${metis_graph_files}
    #mv ${metis_graph_files}/event_graph.graph ${metis_graph_files}/event_graph_2.graph

    # Produce partition file with METIS
    cd ${metis_path}
    ./${metis_bin} ${metis_graph_files}/${graph_base_1}.graph ${n_procs}
    ./${metis_bin} ${metis_graph_files}/${graph_base_2}.graph ${n_procs}
    part_file_1=${graph_base_1}.graph.part.${n_procs}
    part_file_2=${graph_base_2}.graph.part.${n_procs}
    mv ${metis_graph_files}/${part_file_1} ${part_dir}/${part_file_1}
    mv ${metis_graph_files}/${part_file_2} ${part_dir}/${part_file_2}
    cd -

fi
if [ ${partition} == "divide" ]; then

    part_file_1=${graph_base_1}.graph.part.${n_procs}
    part_file_2=${graph_base_2}.graph.part.${n_procs}
    mpirun -np ${n_procs} ${static_assign} ${input_graph1} ${input_graph2} "mod_divide" ${part_dir}/${part_file_1} ${part_dir}/${part_file_2} 	

fi

# Perform Runs of Fido
for run_idx in `seq -f "%03g" ${run_idx_low} ${run_idx_high}`; do
    
        # Prepare output path
	echo ${results_path}
	run_path=${results_path}/sims_$LSB_JOBINDEX/procs_${n_procs}/run_${run_idx}/
    	mkdir -p ${run_path}
    	cd ${run_path}
        mkdir -p ${run_path}/runtime_data

    	echo "======================================================================"
    	echo "Starting run ${run_idx} of ${input_graph1} with ${input_graph2} on ${n_procs} processes"
    	#mpirun -np ${n_procs} > ${run_path}/graph_align_out.txt 2> ${run_path}/graph_align_err.txt valgrind --leak-check=full --error-limit=no --log-file="valgrind_out.txt" --suppressions=${val_mpi_suppr2} --suppressions=${val_mpi_suppr} --suppressions=/home/pnbell/Src_GraphAlignment/GRAAL/parallel_graph_align/mpi_supp_samp.supp --gen-suppressions=all ${g_align} ${input_graph1} ${input_graph2} ${orbit_file} ${time_keeping}
	#mpirun -np ${n_procs} > ${run_path}/graph_align_out.txt 2> ${run_path}/graph_align_err.txt valgrind --leak-check=full --error-limit=no --log-file="valgrind_out.txt" --suppressions=${val_mpi_suppr2} --suppressions=${val_mpi_suppr} --suppressions=/home/pnbell/Src_GraphAlignment/GRAAL/parallel_graph_align/mpi_supp_samp.supp ${g_align} ${input_graph1} ${input_graph2} ${orbit_file} ${time_keeping}
	mpirun -np ${n_procs} > ${run_path}/graph_align_out.txt 2> ${run_path}/graph_align_err.txt ${g_align} ${input_graph1} ${input_graph2} ${orbit_file} ${time_keeping} ${part_dir}/${part_file_1} ${part_dir}/${part_file_2}
	#mpirun -np ${n_procs} > ${run_path}/deg_count_out.txt 2> ${run_path}/deg_count_err.txt ${deg_count} ${input_graph1} ${input_graph2}
	#mpirun -np ${n_procs} > graph_align_out.txt 2> graph_align_err.txt ${g_align} ${input_graph1} ${input_graph2} ${orbit_file} ${time_keeping}
    	echo "======================================================================"

done

# Construct deg distribution across processes
if [ "${load_assignment}" == "static" ]; then
    cd ${run_path}/../
    mpirun -np ${n_procs} > deg_count_out.txt 2> deg_count_err.txt ${deg_count} ${input_graph1} ${input_graph2} ${metis_graph_files}/event_graph_1.graph.part.${n_procs} ${metis_graph_files}/event_graph_2.graph.part.${n_procs}
fi



