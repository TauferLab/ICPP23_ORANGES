#!/usr/bin/env bash
#PBS -l select=32:system=polaris
#PBS -l place=scatter:exclhost
#PBS -l walltime=12:00:00
#PBS -l filesystems=home:eagle:grand
#PBS -q prod 
#PBS -A Veloc
#PBS -o output/oranges_gpu_distributed32_polaris.out
#PBS -e output/oranges_gpu_distributed32_polaris.err

#export KOKKOS_PROFILE_LIBRARY=$HOME/kokkos-tools/profiling/memory-events/kp_memory_events.so
#export KOKKOS_PROFILE_LIBRARY=$HOME/kokkos-tools/profiling/simple-kernel-timer/kp_kernel_timer.so

NNODES=`wc -l < $PBS_NODEFILE`
NUM_NODES_PER_MPI=8
NRANKS_PER_NODE=4
NDEPTH=8
NTHREADS=8

TIME_DATA=oranges_time_data.txt
chunk_sizes=(128)

#dedup_approaches=('--run-full-chkpt' '--run-basic-chkpt' '--run-list-chkpt' '--run-tree-chkpt' '--run-tree-low-offset-ref-chkpt' '--run-tree-low-offset-chkpt' '--run-tree-low-root-ref-chkpt' '--run-tree-low-root-chkpt')
#approaches=('--run-full-chkpt' '--run-basic-chkpt' '--run-list-chkpt' '--run-tree-low-offset-chkpt' '--run-tree-low-root-chkpt')
approaches=('--run-full-chkpt' '--run-list-chkpt' '--run-tree-low-offset-chkpt' '--run-tree-low-root-chkpt')
#approaches=('--run-list-chkpt' '--run-tree-low-offset' '--run-tree-low-root')
#approaches1=('--run-full-chkpt' '--run-list-chkpt' )
#approaches2=('--run-tree-low-offset-chkpt' '--run-tree-low-root-chkpt')
#approaches3=('--run-tree-low-root-chkpt')

. load_modules_polaris.sh
module list
nvidia-smi topo -m

NTOTRANKS=$(( NUM_NODES_PER_MPI * NRANKS_PER_NODE ))
echo "NUM_OF_NODES= ${NNODES} NUM_NODES_PER_MPI= ${NUM_NODES_PER_MPI} TOTAL_NUM_RANKS= ${NTOTRANKS} RANKS_PER_NODE= ${NRANKS_PER_NODE} THREADS_PER_RANK= ${NTHREADS}"

cd ${PBS_O_WORKDIR}

#GRAPH1=test_graph_scale_2_25_25_25_25_A_16777216.edgelist
GRAPH1=delaunay_n24/delaunay_n24.mtx
#GRAPH1=nlpkkt200/nlpkkt200.mtx
GRAPH2=germany_osm.mtx
#INTERVAL=2380000000
#INTERVAL=10000
INTERVAL=10000000
MAX_INTERVAL=10
TIME_INTERVAL=600

cd 32proc
counter=1
for i in $(seq 4); do
  hoststr1=$(sed "${counter}q;d" $PBS_NODEFILE)
  echo "$hoststr1" > "./hostfile${i}"
  counter=$((counter+1))
  hoststr2=$(sed "${counter}q;d" $PBS_NODEFILE)
  echo "$hoststr2" >> "./hostfile${i}"
  counter=$((counter+1))
  hoststr3=$(sed "${counter}q;d" $PBS_NODEFILE)
  echo "$hoststr3" >> "./hostfile${i}"
  counter=$((counter+1))
  hoststr4=$(sed "${counter}q;d" $PBS_NODEFILE)
  echo "$hoststr4" >> "./hostfile${i}"
  counter=$((counter+1))
  hoststr5=$(sed "${counter}q;d" $PBS_NODEFILE)
  echo "$hoststr5" >> "./hostfile${i}"
  counter=$((counter+1))
  hoststr6=$(sed "${counter}q;d" $PBS_NODEFILE)
  echo "$hoststr6" >> "./hostfile${i}"
  counter=$((counter+1))
  hoststr7=$(sed "${counter}q;d" $PBS_NODEFILE)
  echo "$hoststr7" >> "./hostfile${i}"
  counter=$((counter+1))
  hoststr8=$(sed "${counter}q;d" $PBS_NODEFILE)
  echo "$hoststr8" >> "./hostfile${i}"
  counter=$((counter+1))
done

for chunk_size in "${chunk_sizes[@]}"
do
  for i in $(seq 5); do
    counter=1
    for approach in "${approaches[@]}" 
    do
      echo "Launching mpiexec w/"
      cat "hostfile${counter}"
      mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --hostfile "hostfile${counter}" --depth=8 --cpu-bind depth ../set_affinity_gpu_polaris.sh ../fido ../test_graphs/$GRAPH1 ../test_graphs/$GRAPH2 ../orbit_signatures.txt $TIME_DATA $INTERVAL $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$NTHREADS --kokkos-map-device-id-by=mpi_rank &
#      mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --depth=8 --cpu-bind depth ../set_affinity_gpu_polaris.sh ../fido ../test_graphs/$GRAPH1 ../test_graphs/$GRAPH2 ../orbit_signatures.txt $TIME_DATA $INTERVAL $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$NTHREADS --kokkos-map-device-id-by=mpi_rank 
      counter=$((counter+1))
      sleep 1s
    done
    wait
  done
done


#cd 8proc
#for chunk_size in "${chunk_sizes[@]}"
#do
#  counter=1
#  for approach in "${approaches1[@]}" 
#  do
#    for i in $(seq 5); do
#      echo "Launching mpiexec w/"
#      cat "hostfile${i}"
#      mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --hostfile "hostfile{counter}" --depth=8 --cpu-bind depth ../set_affinity_gpu_polaris.sh ../fido ../test_graphs/$GRAPH1 ../test_graphs/$GRAPH2 ../orbit_signatures.txt $TIME_DATA $INTERVAL $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$NTHREADS --kokkos-map-device-id-by=mpi_rank &
#      counter=$((counter+1))
#      sleep 1s
#    done
#  done
#  wait
#  counter=1
#  for approach in "${approaches2[@]}" 
#  do
#    for i in $(seq 5); do
#      echo "Launching mpiexec w/ $hoststr1 and $hoststr2"
#      cat "hostfile${i}"
#      mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --hostfile "hostfile{counter}" --depth=8 --cpu-bind depth ../set_affinity_gpu_polaris.sh ../fido ../test_graphs/$GRAPH1 ../test_graphs/$GRAPH2 ../orbit_signatures.txt $TIME_DATA $INTERVAL $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$NTHREADS --kokkos-map-device-id-by=mpi_rank &
#      counter=$((counter+1))
#      sleep 1s
#    done
#  done
#  wait
#  counter=1
#  for approach in "${approaches3[@]}" 
#  do
#    for i in $(seq 5); do
#      echo "Launching mpiexec w/ $hoststr1 and $hoststr2"
#      cat "hostfile${i}"
#      mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --hostfile "hostfile{counter}" --depth=8 --cpu-bind depth ../set_affinity_gpu_polaris.sh ../fido ../test_graphs/$GRAPH1 ../test_graphs/$GRAPH2 ../orbit_signatures.txt $TIME_DATA $INTERVAL $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$NTHREADS --kokkos-map-device-id-by=mpi_rank &
#      counter=$((counter+1))
#      sleep 1s
#    done
#  done
#  wait
#done



