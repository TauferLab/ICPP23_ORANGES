#!/usr/bin/bash

NUM_PROCESSES=1
CORES_PER_PROC=16
THREADS_PER_PROC=2

GRAPHS=(
	"/home/cc/icpp_graphs/message_race_11174336_Gorder5_updated.edgelist" 
	"/home/cc/icpp_graphs/unstructured_mesh_14418368_Gorder5_updated.edgelist" 
	"/home/cc/icpp_graphs/asia_osm_Gorder5_updated.mtx" 
	"/home/cc/icpp_graphs/hugebubbles-00000_Gorder5_updated.edgelist"
)
INTERVALS=(
	12500000000 
	23000000000 
	810000000 
	5700000000
)

TIME_DATA=time_data.txt

MAX_INTERVAL=5
TIME_INTERVAL=0

num_iter=3

chunksizes=(128 256 512)

approaches=(
	'--run-full-chkpt' 
	'--run-basic-chkpt' 
	'--run-list-chkpt' 
	'--run-tree-low-offset-chkpt'
)

#export OMP_PLACES=threads
#export OMP_PROC_BIND=spread


for chunk_size in "${chunksizes[@]}"
do
  mkdir -p icpp_data/message_race/vary_chunk_size/dedup/$chunk_size
  mkdir -p icpp_data/message_race/vary_chunk_size/dedup/$chunk_size
  mkdir -p icpp_data/message_race/vary_chunk_size/dedup/$chunk_size
  mkdir -p icpp_data/unstructured_mesh/vary_chunk_size/dedup/$chunk_size
  mkdir -p icpp_data/unstructured_mesh/vary_chunk_size/dedup/$chunk_size
  mkdir -p icpp_data/unstructured_mesh/vary_chunk_size/dedup/$chunk_size
  mkdir -p icpp_data/asia_osm/vary_chunk_size/dedup/$chunk_size
  mkdir -p icpp_data/asia_osm/vary_chunk_size/dedup/$chunk_size
  mkdir -p icpp_data/asia_osm/vary_chunk_size/dedup/$chunk_size
  mkdir -p icpp_data/hugebubbles/vary_chunk_size/dedup/$chunk_size
  mkdir -p icpp_data/hugebubbles/vary_chunk_size/dedup/$chunk_size
  mkdir -p icpp_data/hugebubbles/vary_chunk_size/dedup/$chunk_size

  for i in "${!GRAPHS[@]}"
  do
    GRAPH=${GRAPHS[i]}
    echo ${GRAPH}
    INTERVAL=${INTERVALS[i]}
    echo ${INTERVAL}
    IFS='\/' read -ra GRAPH_NAME <<< "$GRAPH"
    logname=${GRAPH_NAME[-1]}
    for approach in "${approaches[@]}"
    do
      for iter in $(seq 1 $num_iter)
      do
        max=0
        echo "Approach: $approach"
        mpirun -n $NUM_PROCESSES $HOME/Src_ORANGES/build/oranges "${GRAPH}" $HOME/Src_ORANGES/data/orbit_signatures.txt $TIME_DATA ${INTERVAL} $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$THREADS_PER_PROC --kokkos-map-device-id-by=mpi_rank &
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
    done
    sed -i " 1 s/.*/&,Max GPU Memory (B)/" ${logname}*${chunk_size}.csv
    if [[ ${logname} == *"message_race"* ]]; then
      mv ${logname}*.csv icpp_data/message_race/vary_chunk_size/dedup/${chunk_size}/
    elif [[ ${logname} == *"unstructured_mesh"* ]]; then
      mv ${logname}*.csv icpp_data/unstructured_mesh/vary_chunk_size/dedup/${chunk_size}/
    elif [[ ${logname} == *"asia_osm"* ]]; then
      mv ${logname}*.csv icpp_data/asia_osm/vary_chunk_size/dedup/${chunk_size}/
    elif [[ ${logname} == *"hugebubbles"* ]]; then
      mv ${logname}*.csv icpp_data/hugebubbles/vary_chunk_size/dedup/${chunk_size}/
    fi
  done
done


