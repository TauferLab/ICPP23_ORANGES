chkpt_dir=$HOME/checkpoint_files/unstructured_mesh_gorder_full_chkpts
graph_name=unstructured_mesh_14418368_Gorder_updated.edgelist

last_id=5

# Cascaded
echo "--------------------"
echo "------Cascaded------"
echo "--------------------"
for i in $(seq 0 $last_id); do
  filename="${graph_name}_Rank0.${i}.hashtree.incr_chkpt.full.chkpt"
  logname="${filename}.cascaded.csv"
  ./bin/benchmark_cascaded_chunked -c true -i 1 -f "${chkpt_dir}/${filename}" > $logname &
  PROC_ID=$!
  max=0
  while kill -0 "$PROC_ID" >/dev/null 2>&1; do
    curr=$(nvidia-smi --query-gpu=memory.used --format=csv|grep -v memory|awk '{print $1}' 2>&1 | head -n 1)
    [ $curr -gt $max ] && max=$curr
    sleep .05
  done
  max=$(($max * 1024 * 1024))
  sed -i " 1 s/.*/&,Max GPU Memory (B)/" $logname
  sed -i " 2 s/.*/&,${max}/" $logname
done

# Deflate
echo "--------------------"
echo "-------Deflate------"
echo "--------------------"
for i in $(seq 0 $last_id); do
  filename="${graph_name}_Rank0.${i}.hashtree.incr_chkpt.full.chkpt"
  logname="${filename}.deflate.csv"
  ./bin/benchmark_deflate_chunked -c true -i 1 -f "${chkpt_dir}/${filename}" > $logname &
  PROC_ID=$!
  max=0
  while kill -0 "$PROC_ID" >/dev/null 2>&1; do
    curr=$(nvidia-smi --query-gpu=memory.used --format=csv|grep -v memory|awk '{print $1}' 2>&1 | head -n 1)
    [ $curr -gt $max ] && max=$curr
    sleep .05
  done
  max=$(($max * 1024 * 1024))
  sed -i " 1 s/.*/&,Max GPU Memory (B)/" $logname
  sed -i " 2 s/.*/&,${max}/" $logname
done

# GDeflate
echo "--------------------"
echo "------GDeflate------"
echo "--------------------"
for i in $(seq 0 $last_id); do
  filename="${graph_name}_Rank0.${i}.hashtree.incr_chkpt.full.chkpt"
  logname="${filename}.gdeflate.csv"
  ./bin/benchmark_gdeflate_chunked -c true -i 1 -f "${chkpt_dir}/${filename}" > $logname &
  PROC_ID=$!
  max=0
  while kill -0 "$PROC_ID" >/dev/null 2>&1; do
    curr=$(nvidia-smi --query-gpu=memory.used --format=csv|grep -v memory|awk '{print $1}' 2>&1 | head -n 1)
    [ $curr -gt $max ] && max=$curr
    sleep .05
  done
  max=$(($max * 1024 * 1024))
  sed -i " 1 s/.*/&,Max GPU Memory (B)/" $logname
  sed -i " 2 s/.*/&,${max}/" $logname
done

# LZ4
echo "--------------------"
echo "---------LZ4--------"
echo "--------------------"
for i in $(seq 0 $last_id); do
  filename="${graph_name}_Rank0.${i}.hashtree.incr_chkpt.full.chkpt"
  logname="${filename}.lz4.csv"
  ./bin/benchmark_lz4_chunked -c true -i 1 -f "${chkpt_dir}/${filename}" > $logname &
  PROC_ID=$!
  max=0
  while kill -0 "$PROC_ID" >/dev/null 2>&1; do
    curr=$(nvidia-smi --query-gpu=memory.used --format=csv|grep -v memory|awk '{print $1}' 2>&1 | head -n 1)
    [ $curr -gt $max ] && max=$curr
    sleep .05
  done
  max=$(($max * 1024 * 1024))
  sed -i " 1 s/.*/&,Max GPU Memory (B)/" $logname
  sed -i " 2 s/.*/&,${max}/" $logname
done

# Snappy
echo "--------------------"
echo "-------Snappy-------"
echo "--------------------"
for i in $(seq 0 $last_id); do
  filename="${graph_name}_Rank0.${i}.hashtree.incr_chkpt.full.chkpt"
  logname="${filename}.snappy.csv"
  ./bin/benchmark_snappy_chunked -c true -i 1 -f "${chkpt_dir}/${filename}" > $logname &
  PROC_ID=$!
  max=0
  while kill -0 "$PROC_ID" >/dev/null 2>&1; do
    curr=$(nvidia-smi --query-gpu=memory.used --format=csv|grep -v memory|awk '{print $1}' 2>&1 | head -n 1)
    [ $curr -gt $max ] && max=$curr
    sleep .05
  done
  max=$(($max * 1024 * 1024))
  sed -i " 1 s/.*/&,Max GPU Memory (B)/" $logname
  sed -i " 2 s/.*/&,${max}/" $logname
done

# Zstd
echo "--------------------"
echo "--------Zstd--------"
echo "--------------------"
for i in $(seq 0 $last_id); do
  filename="${graph_name}_Rank0.${i}.hashtree.incr_chkpt.full.chkpt"
  logname="${filename}.zstd.csv"
  ./bin/benchmark_zstd_chunked -c true -i 1 -f "${chkpt_dir}/${filename}" > $logname &
  PROC_ID=$!
  max=0
  while kill -0 "$PROC_ID" >/dev/null 2>&1; do
    curr=$(nvidia-smi --query-gpu=memory.used --format=csv|grep -v memory|awk '{print $1}' 2>&1 | head -n 1)
    [ $curr -gt $max ] && max=$curr
    sleep .05
  done
  max=$(($max * 1024 * 1024))
  sed -i " 1 s/.*/&,Max GPU Memory (B)/" $logname
  sed -i " 2 s/.*/&,${max}/" $logname
done


