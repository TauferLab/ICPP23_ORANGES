#!/usr/bin/bash

NUM_PROCESSES=1
CORES_PER_PROC=16
THREADS_PER_PROC=2

GRAPH=/home/cc/icpp_graphs/delaunay_n24_Gorder5_updated.edgelist
INTERVAL=10000000
SCALES=(1 2 4 8 16 32 64)

TIME_DATA=time_data.txt

MAX_INTERVAL=10
TIME_INTERVAL=0

num_iter=3

chunksize=128

approaches=(
	'--run-full-chkpt' 
	'--run-tree-low-offset-chkpt'
)

#export OMP_PLACES=threads
#export OMP_PROC_BIND=spread


for scale in "${SCALES[@]}"
do
  mkdir -p icpp_data/delaunay_n24/${scale}proc

  echo ${GRAPH}
  echo ${INTERVAL}
  IFS='\/' read -ra GRAPH_NAME <<< "$GRAPH"
  logname=${GRAPH_NAME[-1]}
  for approach in "${approaches[@]}"
  do
    for iter in $(seq 1 $num_iter)
    do
      max=0
      echo "Approach: $approach"
      mpirun -n $scale $HOME/Src_ORANGES/build/oranges "${GRAPH}" $HOME/Src_ORANGES/data/orbit_signatures.txt $TIME_DATA ${INTERVAL} $chunk_size $MAX_INTERVAL $TIME_INTERVAL $approach --kokkos-num-threads=$THREADS_PER_PROC --kokkos-map-device-id-by=mpi_rank &
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



