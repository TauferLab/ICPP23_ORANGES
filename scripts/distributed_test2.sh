#!/usr/bin/env bash
#PBS -l select=10:system=polaris
#PBS -l place=scatter:exclhost
#PBS -l walltime=6:00:00
#PBS -l filesystems=home:eagle:grand
#PBS -q prod
#PBS -A Veloc
#PBS -o output/oranges_gpu_distributed2_polaris.out
#PBS -e output/oranges_gpu_distributed2_polaris.err

#export KOKKOS_PROFILE_LIBRARY=$HOME/kokkos-tools/profiling/memory-events/kp_memory_events.so
#export KOKKOS_PROFILE_LIBRARY=$HOME/kokkos-tools/profiling/simple-kernel-timer/kp_kernel_timer.so

NNODES=`wc -l < $PBS_NODEFILE`
NUM_NODES_PER_MPI=2
NRANKS_PER_NODE=1
NDEPTH=8
NTHREADS=8

TIME_DATA=oranges_time_data.txt
chunk_sizes=(128)

#dedup_approaches=('--run-full-chkpt' '--run-basic-chkpt' '--run-list-chkpt' '--run-tree-chkpt' '--run-tree-low-offset-ref-chkpt' '--run-tree-low-offset-chkpt' '--run-tree-low-root-ref-chkpt' '--run-tree-low-root-chkpt')
#approaches=('--run-full-chkpt' '--run-basic-chkpt' '--run-list-chkpt' '--run-tree-low-offset-chkpt' '--run-tree-low-root-chkpt')
#approaches=('--run-list-chkpt' '--run-tree-low-offset-ref-chkpt' '--run-tree-low-offset-chkpt' '--run-tree-low-root-ref-chkpt' '--run-tree-low-root-chkpt')
#approaches1=('--run-full-chkpt' '--run-basic-chkpt' )
#approaches2=('--run-list-chkpt' '--run-tree-low-offset-chkpt')
#approaches3=('--run-tree-low-root-chkpt')
#approaches1=('--run-tree-low-offset-chkpt' '--run-tree-low-root-chkpt')
approaches=('--run-full-chkpt' '--run-list-chkpt' '--run-tree-low-offset-chkpt' '--run-tree-low-root-chkpt')

. load_modules_polaris.sh
module list
nvidia-smi topo -m

NTOTRANKS=$(( NUM_NODES_PER_MPI * NRANKS_PER_NODE ))
echo "NUM_OF_NODES= ${NNODES} NUM_NODES_PER_MPI= ${NUM_NODES_PER_MPI} TOTAL_NUM_RANKS= ${NTOTRANKS} RANKS_PER_NODE= ${NRANKS_PER_NODE} THREADS_PER_RANK= ${NTHREADS}"

cd ${PBS_O_WORKDIR}

#GRAPH1=test_graph_scale_2_25_25_25_25_A_16777216.edgelist
#GRAPH1=nlpkkt200/nlpkkt200.mtx
GRAPH1=delaunay_n24/delaunay_n24.mtx
GRAPH2=germany_osm.mtx
#INTERVAL=2380000000
INTERVAL=10000000
#INTERVAL=10000
MAX_INTERVAL=10
TIME_INTERVAL=600

#hoststr1=$(sed "1q;d" $PBS_NODEFILE)

#cd 2proc
#for chunk_size in "${chunk_sizes[@]}"
#do
#  for i in $(seq 5)
#  do
#    counter=1
#    for approach in "${approaches[@]}" 
#    do
#      mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --hostfile "hostfile${counter}" --depth=8 --cpu-bind depth ../set_affinity_gpu_polaris.sh ../fido ../test_graphs/$GRAPH1 ../test_graphs/$GRAPH2 ../orbit_signatures.txt $TIME_DATA $INTERVAL $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$NTHREADS --kokkos-map-device-id-by=mpi_rank &
#      sleep 1s
#      counter=$((counter+1))
#    done
#    wait
#  done
#done

# 1 rank per node
cd 2proc_sep
echo $(sed "1q;d" $PBS_NODEFILE) >  "./hostfile1"
echo $(sed "2q;d" $PBS_NODEFILE) >> "./hostfile1"
echo $(sed "3q;d" $PBS_NODEFILE) >  "./hostfile2"
echo $(sed "4q;d" $PBS_NODEFILE) >> "./hostfile2"
echo $(sed "5q;d" $PBS_NODEFILE) >  "./hostfile3"
echo $(sed "6q;d" $PBS_NODEFILE) >> "./hostfile3"
echo $(sed "7q;d" $PBS_NODEFILE) >  "./hostfile4"
echo $(sed "8q;d" $PBS_NODEFILE) >> "./hostfile4"
echo $(sed "9q;d" $PBS_NODEFILE) >  "./hostfile5"
echo $(sed "10q;d" $PBS_NODEFILE) >> "./hostfile5"
for chunk_size in "${chunk_sizes[@]}"
do
  for approach in "${approaches[@]}" 
  do
    for i in $(seq 5); do
      echo "Launching mpiexec w/ $hoststr{i}"
      mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --hostfile "hostfile${i}" --depth=8 --cpu-bind depth ../set_affinity_gpu_polaris.sh ../fido ../test_graphs/$GRAPH1 ../test_graphs/$GRAPH2 ../orbit_signatures.txt $TIME_DATA $INTERVAL $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$NTHREADS --kokkos-map-device-id-by=mpi_rank &
      counter=$((counter+1))
      sleep 1s
    done
    wait
  done

## Multiple ranks per node
#cd 2proc
#counter=1
#for chunk_size in "${chunk_sizes[@]}"
#do
#  for approach in "${approaches1[@]}" 
#  do
#    for i in $(seq 5); do
#      hoststr=$(sed "${counter}q;d" $PBS_NODEFILE)
#      echo "$hoststr" > "./hostfile${counter}"
#      echo "Launching mpiexec w/ $hoststr"
#      mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --hostfile "hostfile${counter}" --depth=8 --cpu-bind depth ../set_affinity_gpu_polaris.sh ../fido ../test_graphs/$GRAPH1 ../test_graphs/$GRAPH2 ../orbit_signatures.txt $TIME_DATA $INTERVAL $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$NTHREADS --kokkos-map-device-id-by=mpi_rank &
#      counter=$((counter+1))
#      sleep 1s
#    done
#    wait
#  done
#  counter=1
#  for approach in "${approaches2[@]}" 
#  do
#    for i in $(seq 5); do
#      hoststr=$(sed "${counter}q;d" $PBS_NODEFILE)
#      echo "$hoststr" > "./hostfile${counter}"
#      echo "Launching mpiexec w/ $hoststr"
#      mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --hostfile "hostfile${counter}" --depth=8 --cpu-bind depth ../set_affinity_gpu_polaris.sh ../fido ../test_graphs/$GRAPH1 ../test_graphs/$GRAPH2 ../orbit_signatures.txt $TIME_DATA $INTERVAL $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$NTHREADS --kokkos-map-device-id-by=mpi_rank &
#      counter=$((counter+1))
#      sleep 1s
#    done
#    wait
#  done
#  counter=1
#  for approach in "${approaches3[@]}" 
#  do
#    for i in $(seq 5); do
#      hoststr=$(sed "${counter}q;d" $PBS_NODEFILE)
#      echo "$hoststr" > "./hostfile${counter}"
#      echo "Launching mpiexec w/ $hoststr"
#      mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --hostfile "hostfile${counter}" --depth=8 --cpu-bind depth ../set_affinity_gpu_polaris.sh ../fido ../test_graphs/$GRAPH1 ../test_graphs/$GRAPH2 ../orbit_signatures.txt $TIME_DATA $INTERVAL $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$NTHREADS --kokkos-map-device-id-by=mpi_rank &
#      counter=$((counter+1))
#      sleep 1s
#    done
#    wait
#  done
done


