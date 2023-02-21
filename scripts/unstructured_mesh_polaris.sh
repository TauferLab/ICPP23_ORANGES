#!/usr/bin/env bash
#PBS -l select=10:system=polaris
#PBS -l place=scatter:exclhost
#PBS -l walltime=5:00:00
#PBS -l filesystems=home:eagle:grand
#PBS -q prod
#PBS -A Veloc
#PBS -o output/oranges_gpu_unstructured_mesh_polaris.out
#PBS -e output/oranges_gpu_unstructured_mesh_polaris.err

#export KOKKOS_PROFILE_LIBRARY=$HOME/kokkos-tools/profiling/memory-events/kp_memory_events.so
#export KOKKOS_PROFILE_LIBRARY=$HOME/kokkos-tools/profiling/simple-kernel-timer/kp_kernel_timer.so

NNODES=`wc -l < $PBS_NODEFILE`
NUM_NODES_PER_MPI=1
NRANKS_PER_NODE=1
NDEPTH=32
NTHREADS=8

TIME_DATA=oranges_time_data.txt
chunk_sizes=(128 256 512 1024 2048 4096)

#dedup_approaches=('--run-full-chkpt' '--run-basic-chkpt' '--run-list-chkpt' '--run-tree-chkpt' '--run-tree-low-offset-ref-chkpt' '--run-tree-low-offset-chkpt' '--run-tree-low-root-ref-chkpt' '--run-tree-low-root-chkpt')
approaches=('--run-full-chkpt' '--run-basic-chkpt' '--run-list-chkpt' '--run-tree-low-offset-chkpt' '--run-tree-low-root-chkpt')
num_iter=1

. load_modules_polaris.sh
module list
nvidia-smi

NTOTRANKS=$(( NUM_NODES_PER_MPI * NRANKS_PER_NODE ))
echo "NUM_OF_NODES= ${NNODES} NUM_NODES_PER_MPI= ${NUM_NODES_PER_MPI} TOTAL_NUM_RANKS= ${NTOTRANKS} RANKS_PER_NODE= ${NRANKS_PER_NODE} THREADS_PER_RANK= ${NTHREADS}"

cd ${PBS_O_WORKDIR}

GRAPH1=unstructured_mesh_nprocs_32_packed_nd_fraction_0.25_niters_2048_msg_size_32_run_000_nodes_14418368.edgelist
GRAPH2=unstructured_mesh_nprocs_32_packed_nd_fraction_0.25_niters_2048_msg_size_32_run_001_nodes_14418368.edgelist
INTERVAL=12000000000 # Unstructured mesh 14418368
MAX_INTERVAL=10
TIME_INTERVAL=0
for chunk_size in "${chunk_sizes[@]}"
do
  for approach in "${approaches[@]}" 
  do
    for i in $(seq "5"); do
      hoststr=$(sed "${i}q;d" $PBS_NODEFILE)
      echo "$hoststr" > "./hostfile${i}"
      echo "Launching mpiexec w/ $hoststr"
      mpiexec -n $NUM_NODES_PER_MPI --ppn ${NRANKS_PER_NODE} --hostfile "hostfile${i}" ./fido ./test_graphs/$GRAPH1 ./test_graphs/$GRAPH2 ./orbit_signatures.txt $TIME_DATA $INTERVAL $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$NTHREADS --kokkos-map-device-id-by=random &
      sleep 1s
    done
#    wait
    for i in $(seq 6 10); do
      hoststr=$(sed "${i}q;d" $PBS_NODEFILE)
      echo "$hoststr" > "./hostfile${i}"
      echo "Launching mpiexec w/ $hoststr"
      mpiexec -n $NUM_NODES_PER_MPI --ppn ${NRANKS_PER_NODE} --hostfile "hostfile${i}" ./fido ./test_graphs/$GRAPH2 ./test_graphs/$GRAPH1 ./orbit_signatures.txt $TIME_DATA $INTERVAL $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$NTHREADS --kokkos-map-device-id-by=random &
      sleep 1s
    done
    wait
  done
done



