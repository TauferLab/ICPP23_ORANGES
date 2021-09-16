#!/usr/bin/env bash


Help() {

	echo ""
	echo "The following command line switches can be used to define parameters for your job submission:"
	echo ""
	echo "[-p]    Defines the number of mpi processes used when parallelizing graph alignment. (Default 2 MPI processes)"
	echo "[-n]    The number of compute nodes requested for parallelizing graph alignment. (Default 1 node)"
	echo "        If you're running on an unscheduled system, this value should be set to 1."
	echo "[-sc]   Used to define which schedule system is currently in use. (Code will request this if not set)"
	echo "        Must be one of the following options: lsf, slurm, or none."
	echo "[-g]    Takes two inputs as arguments corresponding to the full paths to 2 graph files in the edge list format."
	echo "        (Defaults to sgraph1.txt and sgraph2.txt within the test_graphs directory of Fido."
	echo "[-o]    If used, allows the user to define their own path to store output from the project. (Defaults to the directory '$HOME/fido_output/graph_align_<current date and time>/')"
	echo "        When using this flag, be sure to provide an absolute path that can exist on your machine."
	echo "        If you run this script multiple times on the same settings, be sure to use different paths to avoid overlap and overwriting of files."
	echo "[-q]    Defines the queue to submit jobs to when on a scheduled system. (i.e. the Slurm queue or the LSF queue)"
	echo "        (Defaults to the 'normal' queue)"
	echo "[-t]    Defines the maximum time limit set for the job of graph alignment. (Default 10 minutes)"
	echo "[-h]    Used to display this list of switch options."
	echo ""

}


while [ -n "$1" ]; do
	case "$1" in
		-p)  n_procs=$2; shift; shift ;;
		-n)  n_nodes=$2; shift; shift ;;
		-q)  queue=$2; shift; shift ;;
		-t)  time_limit=$2; shift; shift;;
		-g)  graph_one=$2; graph_two=$3; shift; shift; shift ;;
		-o)  results_path=$2; shift; shift ;;
		-sc) scheduler=$2; shift; shift ;;
		-h)  Help; exit ;;
		*)   echo "$1 is not an option"; exit ;;
	esac
done


# Pick a scheduler
while true; do
    #read -p "Which job scheduler would you like to use? Input is case sensitive. (lsf, slurm, unscheduled) " scheduler
	case ${scheduler} in
		"lsf" | "slurm" | "none" ) break ;;
		* ) echo "The scheduler system needs to be selected from one of the three optional schedulers. (lsf, slurm, none) "
		read -p "Please respond with one of the listed options. Input is case sensitive. " scheduler ;;
	esac
done


# Assign needed paths
paths_dir=$(pwd)/fido/parallel_graph_align/
source ${paths_dir}/fido_paths.config
fido_sched_script=${fido_project_root}/fido/parallel_graph_align/graph_align_sched.sh

# Assign Default Values
n_procs="${n_procs:=2}"
scheduler="${scheduler:="none"}"
n_nodes="${n_nodes:=1}"
graph_one="${graph_one:="$(pwd)/fido/test_graphs/sgraph1.txt"}"
graph_two="${graph_two:="$(pwd)/fido/test_graphs/sgraph2.txt"}"
results_path="${results_path:=$HOME/fido_output/graph_align_$(date +%s)/}"
queue="${queue:="normal"}"
time_limit="${time_limit:=10}"
load_balancer="${load_balancer:="dynamic"}"


# Ensure the input values are valid to program requirements
while (( $(echo "$n_procs < 1" |bc -l) )) || ! [[ "$n_procs" =~ ^[0-9]+$ ]] || [ -z "$n_procs" ] ; do
	echo "Number of MPI processes was set too low or is not an integer."
	echo "Please set number of processes to an integer greater than 0."
	read -p "Number of MPI processes requested: " n_procs
done
while (( $(echo "$n_nodes < 1" |bc -l) )) || ! [[ "$n_nodes" =~ ^[0-9]+$ ]] || [ -z "$n_nodes" ] ; do
	echo "Number of compute nodes was set too low or is not an integer."
	echo "Please set number of compute nodes to an integer greater than 0. We recommend using at least 2 if available."
	read -p "Number of compute nodes requested: " n_nodes
done

# Ensure the graph files are valid full paths with files.
while true; do
	case "${graph_one}" in
		/*) 	if [ -f "${graph_one}" ] ; then
				break;
			else
				echo "The resquested input graph 1: ${graph_one} is not a file."
				read -p "Please input a valid filename: " graph_one
			fi ;;
		*)	echo "The requested input graph 1: ${graph_one} is not an absolute path. Please make sure to input an absolute path for your graph files. It should start with a '/' character."
			echo "As an example, you might input something of the form "$(pwd)"/example_graph_file"
			read -p "Input your file here: " graph_one ;;
	esac
done
while true; do
        case "${graph_two}" in
                /*)     if [ -f "${graph_two}" ] ; then
                                break;
                        else
                                echo "The resquested input graph 2: ${graph_two} is not a file."
                                read -p "Please input a valid filename: " graph_two
                        fi ;;
                *)      echo "The requested input graph 2: ${graph_two} is not an absolute path. Please make sure to input an absolute path for your graph files. It should start with a '/' character."
                        echo "As an example, you might input something of the form "$(pwd)"/example_graph_file"
                        read -p "Input your file here: " graph_two ;;
	esac
done

# Ensure the output path is a valid path.
project_path=$(pwd)
cd
while true; do
	case "${results_path}" in
		/*) mkdir -p ${results_path} 
			if [[ -d "${results_path}" ]]; then
				break
			else
				echo "Your output path is not a valid path.  Please make sure that the path you provided can exist on your machine."
				echo "As an example, you might input something of the form "$(pwd)"/example_path"
				read -p "Input your path here: " results_path
			fi ;;
		*)  echo "Your output path is not an absolute path.  Please make sure to input an absolute path for your output. It should start with a '/' character."
			echo "As an example, you might input something of the form "$(pwd)"/example_path" 
			read -p "Input your path here: " results_path ;;
	esac
done
cd ${project_path}


## Start fido jobs
bash ${fido_sched_script} ${n_procs} ${n_nodes} ${graph_one} ${graph_two} ${paths_dir} ${results_path} ${scheduler} ${queue} ${time_limit} ${load_balancer}


echo ""
echo "==========================================================="
echo ""
echo "Your output will be stored in the directory: ${results_path}"
echo ""
echo "==========================================================="
echo ""
