# Fido

## **Software Overview**

This repository contains a software tool for pairwise, one-to-one alignment of graphs. This procedure is broken down into the following two stages:
1. Node Pair Similarity Scoring
2. Node Alignment

### Stage 1 - Node Pair Similarity Scoring

Stage 1 takes as input 2 graphs in the '.graphml' format to align and evaluate the similarity of each pair of nodes across those graphs.  This similarity is constructed using the following steps:
1. Calculate the graphlet degree vector (GDV) for each node of both input graphs.
2. Calculate the L2 norm similarity between the GDVs for each node across the input the input graphs.

This similarity matrix is stored for use in the second stage, node alignment, of Fido.

Code for stage 1 is found within directories `fido/parallel_graph_align/main_code_parallel.cpp` and `fido/headers/`.

### Stage 2 - Node Alignment

Stage 2 is currently yet to be implemented.

## **Installation**

Use the following steps to prepare and install Fido:
1. Make sure that a **C compiler** and a **version of MPI** are installed on the system you're running with prior to starting any Fido specific installations or submission of scripts.
2. Once you have a **C compiler** and a **version of MPI** installed, clone the git repository of Fido onto your local machine.
3. Enter the projects root directory and use the following command to install the rest of Fidos dependencies along with Fido itself:

```
. ./setup_fido.sh
```

How to setup Fido with oneAPI:
1. Go to Download the Intel oneAPI Base Toolkit to retrieve the oneAPI base toolkit. Select the correct OS and package manager for your system and click on the link titled “View Installation Instructions”
2. Follow the guide that is written out. After installing the oneAPI Base Toolkit, install the Intel oneAPI HPC Toolkit. Instructions for this should be on the same page after scrolling down.
3. To compile with DPCPP, add the following line to the top of the Src_Fido/CMakeLists.txt and the Src_Fido/fido/CmakeLists.txt files
* set(CMAKE_CXX_COMPILER “dpcpp”)

## Dependencies

Below is a list of dependencies for running Fido.  If you follow the steps in the section above on "Installation", each of these should be set up on your system.

### Installed by User:
The following packages must be installed prior to running the script 'setup_fido.sh' during project installation.
* A C Compiler (e.g. GCC, ICC)
* A version of MPI (e.g. openmpi, mvapich2)
* Intel oneAPI Base Toolkit
* Intel oneAPI HPC Toolkit

### Submodule Packages (Installed during setup):
The following packages will be installed automatically via the script 'setup_fido.sh' during project installation.
* [Kokkos](https://github.com/kokkos/kokkos)
* [Kokkos-Kernels](https://github.com/kokkos/kokkos-kernels)

## **Running Fido**

Use the 'graph_align_fido.sh' script to construct similarity matrices for two input graphs.

The following command line switches can be used to define parameters for your job submission:
* -p        : Defines the size of the mpi communicator (number of MPI processes)
                used when parallelizing GDV calculation. 
                (Default 2 MPI processes)
* -n        : The number of compute nodes requested for running Fido. 
                If you're running on an unscheduled system,
                this value should be set to 1.
                (Default 1 node)
* -sc      : Used to define which schedule system is currently in use.
                Must be one of the following options: lsf, slurm, or unscheduled.
* -g       : Takes two inputs as arguments corresponding to the full paths to 
                2 graph files in the edge list format.
* -o        : If used, allows the user to define their own path to 
                store output from the project. 
                Be sure to define an absolute path that can exist on your machine.
                Use a different path when running multiple times on the same settings to avoid overwriting.
                (Defaults to the directory '$HOME/fido_output/graph_align_\<current date and time\>')
* -q        : Defines the queue to submit jobs to when on a scheduled system.
                Specifically, choose a queue when using the Slurm or LSF scheduler.
                (Defaults to the 'normal' queue)
                **Note**, this will not be used when running on an unscheduled system.
* -t        : Defines the maximum time limit set for the job of graph alignment. 
                (Default 10 minutes)
                **Note**, this will not be used when running on an unscheduled system.
* -h        : Used to display this list of switch options.

## **Output of Fido**

You will find all output for Fido within the output directory defined by the command line switch '-o' above.  **Note** If you did not define your own output for the project, it will use the '$HOME/fido_output/graph_align_\<current date and time\>'.

The output of Fido will take the following directory structure:
* Standard output and error files for the scheduler during runtime.
* A directory defined by the settings of your project. (e.g. number of processes, number of compute nodes)
  * A file containing the calculated GDV for each vertex of input graph 1. (out_gdv_1_\<the title of your graph 1 file\>.txt)
  * A file containing the calculated GDV for each vertex of input graph 2. (out_gdv_2_\<the title of your graph 2 file\>.txt)
  * A file containing the matrix of all-to-all similarity scores across both input graphs. (out_similarity_matrix.txt)
  * A directory storing files associated with debugging.
    * Standard output and standard error files for the alignment code.




