#!/usr/bin/env bash
#BSUB -n 2
#BSUB -gpu "num=2:mode=exclusive_process:mps=no"
#BSUB -R "span[ptile=2]"
#BSUB -o output/%J.out
#BSUB -e output/%J.err

NUM_PROCESSES=2
CORES_PER_PROC=1
THREADS_PER_PROC=4

#export KOKKOS_TOOLS_LIBS=/home/ntan1/kokkos-tools/profiling/space-time-stack/kp_space_time_stack.so

#GRAPH1=test_graphs/sgraph1.txt
#GRAPH1=test_graphs/ecology1.mtx
#GRAPH1=test_graphs/asia_osm.mtx
#GRAPH1=test_graphs/message_race_nprocs_32_niters_2048_msg_size_32_run_000_size_1397184.edgelist
#GRAPH1=test_graphs/message_race_nprocs_32_niters_4096_msg_size_32_run_000_size_2793920.edgelist
#GRAPH1=test_graphs/unstructured_mesh_nprocs_32_packed_nd_fraction_0.25_niters_256_msg_size_32_run_000_nodes_1802687.edgelist
#GRAPH1=/data/gclab/fido/message_race_nprocs_32_niters_2048_msg_size_32_run_000_size_1397184.edgelist
GRAPH1=/data/gclab/fido/message_race_nprocs_32_niters_16384_msg_size_32_run_000_size_11174336.edgelist
#GRAPH1=/data/gclab/fido/unstructured_mesh_nprocs_32_packed_nd_fraction_0.5_niters_256_msg_size_32_run_000_nodes_1802687.edgelist
#GRAPH1=/data/gclab/fido/unstructured_mesh_nprocs_32_packed_nd_fraction_0.5_niters_2048_msg_size_32_run_000_nodes_14418368.edgelist
#GRAPH1=$HOME/RMAT_Generator/test_graph_25_25_25_25_A_131072.edgelist
#GRAPH1=$HOME/RMAT_Generator/test_graph_25_25_25_25_A_1048576.edgelist
#GRAPH1=$HOME/RMAT_Generator/test_graph_scale_4_25_25_25_25_A_1048576.edgelist
#GRAPH1=$HOME/RMAT_Generator/test_graph_25_25_25_25_A_16777216.edgelist

#GRAPH2=test_graphs/sgraph1.txt
#GRAPH2=test_graphs/ecology2.mtx
#GRAPH2=test_graphs/germany_osm.mtx
#GRAPH2=test_graphs/message_race_nprocs_32_niters_2048_msg_size_32_run_001_size_1397184.edgelist
#GRAPH2=test_graphs/message_race_nprocs_32_niters_4096_msg_size_32_run_001_size_2793920.edgelist
#GRAPH2=test_graphs/unstructured_mesh_nprocs_32_packed_nd_fraction_0.25_niters_256_msg_size_32_run_001_nodes_1802687.edgelist
#GRAPH2=/data/gclab/fido/message_race_nprocs_32_niters_2048_msg_size_32_run_001_size_1397184.edgelist
GRAPH2=/data/gclab/fido/message_race_nprocs_32_niters_16384_msg_size_32_run_001_size_11174336.edgelist
#GRAPH2=/data/gclab/fido/unstructured_mesh_nprocs_32_packed_nd_fraction_0.5_niters_256_msg_size_32_run_001_nodes_1802687.edgelist
#GRAPH2=/data/gclab/fido/unstructured_mesh_nprocs_32_packed_nd_fraction_0.5_niters_2048_msg_size_32_run_001_nodes_14418368.edgelist
#GRAPH2=$HOME/RMAT_Generator/test_graph_25_25_25_25_B_131072.edgelist
#GRAPH2=$HOME/RMAT_Generator/test_graph_25_25_25_25_B_1048576.edgelist
#GRAPH2=$HOME/RMAT_Generator/test_graph_scale_4_25_25_25_25_B_1048576.edgelist
#GRAPH2=$HOME/RMAT_Generator/test_graph_25_25_25_25_B_16777216.edgelist

TIME_DATA=time_data.txt

#INTERVAL=1000
#INTERVAL=800000000 # Message Race 11174336
#INTERVAL=12000000000 # Unstructured mesh 14418368
INTERVAL=7244620443 # Message race 11174336, 1 rank
#INTERVAL=3622310222 # Message race 11174336, 2 ranks
#INTERVAL=1811155111 # Message race 11174336, 4 ranks
#INTERVAL=905577556 # Message race 11174336, 8 ranks
#INTERVAL=500000000 # Open Street Maps 
#INTERVAL=10000000 # RMAT 1048576
#INTERVAL=12000000000 # RMAT 16777216
CHUNK_SIZE=4096

export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#mpirun -n $NUM_PROCESSES --map-by slot:PE=$CORES_PER_PROC ./fido $GRAPH1 $GRAPH2 ./../data/orbit_signatures.txt $TIME_DATA $INTERVAL $CHUNK_SIZE --kokkos-num-threads=$THREADS_PER_PROC --kokkos-num-devices=2

mpirun -n $NUM_PROCESSES --map-by slot:PE=$CORES_PER_PROC ./fido $GRAPH1 $GRAPH2 ./../data/orbit_signatures.txt $TIME_DATA $INTERVAL $CHUNK_SIZE --kokkos-num-threads=$THREADS_PER_PROC --kokkos-map-device-id-by=mpi_rank


