#!/bin/bash

#choice="list"
#choice="dir"

#file_list=()
project_root=/home/pnbell/Src_GraphAlignment
#file_dir=${project_root}/backup_slices/slices
file_dir=/home/pnbell/Src_ANACIN-X/strong_scaling_tests

cd ${project_root}/GRAAL/strong_scaling

# Convert each graphml to txt in the file director
for file in ${file_dir}/*.graphml; do
    echo ${file}
    python3 ${project_root}/Pre-Processing/Event_graph_dataset_preprocessing.py ${file}
done
