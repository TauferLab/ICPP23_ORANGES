#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -C gpu
#SBATCH -p exclusive
#SBATCH -t 12:00:00
#SBATCH -o output/%J.out
#SBATCH -e output/%J.err

NUM_PROCESSES=1
CORES_PER_PROC=16
THREADS_PER_PROC=4

#export KOKKOS_TOOLS_LIBS=/home/ntan1/kokkos-tools/kp_memory_events.so
#export KOKKOS_TOOLS_LIBS=/home/ntan1/kokkos-tools/profiling/space-time-stack/kp_space_time_stack.so
#export KOKKOS_TOOLS_LIBS=/home/ntan1/kokkos-tools/profiling/simple-kernel-timer/kp_kernel_timer.so

#GRAPH1=test_graphs/sgraph1.txt
#GRAPH1=test_graphs/ecology1.mtx
#GRAPH1=test_graphs/asia_osm.mtx
GRAPH1=/data/gclab/fido/message_race_nprocs_32_niters_16384_msg_size_32_run_000_size_11174336.edgelist
#GRAPH1=/data/gclab/fido/message_race_nprocs_32_niters_32768_msg_size_4096_run_000_size_22348224.edgelist
#GRAPH1=/data/gclab/fido/unstructured_mesh_nprocs_32_packed_nd_fraction_0.5_niters_2048_msg_size_32_run_000_nodes_14418368.edgelist
#GRAPH1=$HOME/RMAT_Generator/test_graph_25_25_25_25_A_131072.edgelist
#GRAPH1=$HOME/RMAT_Generator/test_graph_25_25_25_25_A_1048576.edgelist
#GRAPH1=$HOME/RMAT_Generator/test_graph_scale_4_25_25_25_25_A_131072.edgelist
#GRAPH1=$HOME/RMAT_Generator/test_graph_scale_4_25_25_25_25_A_1048576.edgelist
#GRAPH1=$HOME/RMAT_Generator/test_graph_25_25_25_25_A_16777216.edgelist
#GRAPH1=$HOME/RMAT_Generator/test_graph_scale_2_25_25_25_25_A_16777216.edgelist
#GRAPH1=test_graphs/wikipedia-20060925.mtx
#GRAPH1=test_graphs/nlpkkt200/nlpkkt200.mtx # Too expensive to calc num combinations
#GRAPH1=test_graphs/delaunay_n23/delaunay_n23.mtx
#GRAPH1=test_graphs/delaunay_n24/delaunay_n24.mtx
#GRAPH1=test_graphs/hugebubbles-00000/hugebubbles-00000.mtx
#GRAPH1=test_graphs/hugebubbles-00010/hugebubbles-00010.mtx
#GRAPH1=test_graphs/mawi_201512012345/mawi_201512012345.mtx # Too big
#GRAPH1=test_graphs/road_central/road_central.mtx
#GRAPH1=test_graphs/rgg_n_2_24_s0/rgg_n_2_24_s0.mtx
#GRAPH1=test_graphs/rgg_n_2_23_s0/rgg_n_2_23_s0.mtx
#GRAPH1=test_graphs/rgg_n_2_22_s0/rgg_n_2_22_s0.mtx # Works, very sparse
#GRAPH1=test_graphs/great-britain_osm/great-britain_osm.mtx
#GRAPH1=test_graphs/rel9/rel9.mtx # Too high degree
#GRAPH1=test_graphs/channel-500x100x100-b050/channel-500x100x100-b050.mtx
#GRAPH1=test_graphs/rajat31/rajat31.mtx
#GRAPH1=test_graphs/CL-10M-1d8-L5.edges

GRAPH2=test_graphs/sgraph1.txt

TIME_DATA=time_data.txt

MAX_INTERVAL=5
TIME_INTERVAL=0

#INTERVAL=100
#INTERVAL=10000000
#INTERVAL=100000000 # RGG
#INTERVAL=1000000000 # mawi
#INTERVAL=1000000000 # nlpkkt200
#INTERVAL=4760000000 # Hugebubbles 00000
#INTERVAL=2380000000 # Delaunay N24
INTERVAL=8000000000 # Message Race 11174336
#INTERVAL=40000000000 # Message Race 22348224
#INTERVAL=12000000000 # Unstructured mesh 14418368
#INTERVAL=7244620443 # Message race 11174336, 1 rank
#INTERVAL=3622310222 # Message race 11174336, 2 ranks
#INTERVAL=1811155111 # Message race 11174336, 4 ranks
#INTERVAL=905577556 # Message race 11174336, 8 ranks
#INTERVAL=500000000 # Open Street Maps 
#INTERVAL=10000000 # RMAT 1048576
#INTERVAL=12000000000 # RMAT 16777216
#chunksizes=(128 256 512 1024 2048 4096)
chunksizes=(256)

#dedup_approaches=('--run-full-chkpt' '--run-basic-chkpt' '--run-list-chkpt' '--run-tree-chkpt' '--run-tree-low-offset-ref-chkpt' '--run-tree-low-offset-chkpt' '--run-tree-low-root-ref-chkpt' '--run-tree-low-root-chkpt')
approaches=('--run-full-chkpt' '--run-basic-chkpt' '--run-list-chkpt' '--run-tree-low-offset-chkpt' '--run-tree-low-root-chkpt')

export OMP_PLACES=threads
export OMP_PROC_BIND=spread

for chunk_size in "${chunksizes[@]}"
do
  IFS='\/' read -ra GRAPH_NAME <<< "$GRAPH1"
  logname=${GRAPH_NAME[-1]}
  for approach in "${approaches[@]}" 
  do
    max=0
    echo "Approach: $approach"
    srun -n $NUM_PROCESSES -C gpu ./fido $GRAPH1 $GRAPH2 ./../data/orbit_signatures.txt $TIME_DATA $INTERVAL $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$THREADS_PER_PROC &
    PROC_ID=$!
    while kill -0 "$PROC_ID" >/dev/null 2>&1; do
      curr=$(nvidia-smi --query-gpu=memory.used --format=csv|grep -v memory|awk '{print $1}' 2>&1 | head -n 1)
      [ $curr -gt $max ] && max=$curr
      sleep .001
    done
    max=$(($max * 1024 * 1024))
    echo "Max GPU Mem usage: $max B"
    for FILE in ${logname}*${chunk_size}.csv; do
      if [ -f "$FILE" ]; then
        echo "$(cat $FILE)",$max"" > $FILE
      fi
    done
  done
  sed -i " 1 s/.*/&,Max GPU Memory (B)/" ${logname}*${chunk_size}.csv
done
