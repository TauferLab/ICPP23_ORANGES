#Change these line to a different location if needed
include ./fido_paths.config
proj_root=${fido_project_root}
#root_path=${proj_root}/fido
ESSENS=${proj_root}/ESSENS-master
#ESSENS= ./ESSENS

#Use this to run on Linux
#CXX=g++
#CXXFLAGS=-g -fopenmp 

#Changed CXX flags to run on MAC
#CXX = /usr/local/opt/llvm/bin/clang++
#CXXFLAGS = -I/usr/local/opt/llvm/include -fopenmp
# LDFLAGS = -L/usr/local/opt/llvm/lib

HEADER_PATHS=$(root_path)/headers

CXX=mpicxx
CXXFLAGS=-fopenmp -g -O3

#Including BASIC_IO
BASIC_IO= $(ESSENS)/Core/Basic_IO

BASIC_IO_PREPROCESS=$(BASIC_IO)/PreProcess

BASIC_IO_FORMAT=$(BASIC_IO)/Format
BASIC_IO_FORMAT_INPUT=$(BASIC_IO_FORMAT)/Input
BASIC_IO_FORMAT_OUTPUT=$(BASIC_IO_FORMAT)/Output

BASIC_IO_INC=-I $(BASIC_IO_PREPROCESS)/Exceptions\
 -I $(BASIC_IO_PREPROCESS)/Level0 \
 -I $(BASIC_IO_PREPROCESS)/Level1 \
 -I $(BASIC_IO_PREPROCESS)/Level2 \
 -I $(BASIC_IO_FORMAT_INPUT_FORMAT)/Exceptions \
 -I $(BASIC_IO_FORMAT)/Level0/ \
 -I $(BASIC_IO_FORMAT_INPUT)/Level1/ -I $(BASIC_IO_FORMAT_INPUT)/Level2 \
-I $(BASIC_IO_FORMAT_OUTPUT)/Level1/ -I $(BASIC_IO_FORMAT_OUTPUT)/Level2 


BASIC_SETOP=$(ESSENS)/Core/Basic_SetOps

BASIC_SETOP_INC=-I $(BASIC_SETOP)/Level0 \
-I $(BASIC_SETOP)/Level1/ \
-I $(BASIC_SETOP)/Level2

BASIC_TRAVERSAL=$(ESSENS)/Core/Basic_Traversal

BASIC_TRAVERSAL_INC=-I$(BASIC_TRAVERSAL)/Level0 \
-I$(BASIC_TRAVERSAL)/Level1 \
-I$(BASIC_TRAVERSAL)/Level2

BASIC_CHANGE=$(ESSENS)/Basic_Change

BASIC_CHANGE_INC=-I$(BASIC_CHANGE)/Level0 \
-I$(BASIC_CHANGE)/Level1 \
-I$(BASIC_CHANGE)/Level2

BASIC_ANALYSIS=$(ESSENS)/Basic_Analysis
BASIC_ANALYSIS_INC=-I$(BASIC_ANALYSIS)/Level0

ESSENS_PATHS=$(BASIC_IO_INC) $(BASIC_SETOP_INC) $(BASIC_TRAVERSAL_INC) $(BASIC_CHANGE_INC) $(BASIC_ANALYSIS_INC)
#HEADER_PATHS=$(root_path)/headers

all:    main_code_parallel.cpp $(HEADER_PATHS)/GDV_functions.hpp $(HEADER_PATHS)/class_definitions.hpp $(HEADER_PATHS)/bfs.hpp $(HEADER_PATHS)/print_disconnected_graph.hpp
	echo $(ESSENS_PATHS)
	$(CXX) $(CXXFLAGS) -o graph_alignment $(ESSENS_PATHS) main_code_parallel.cpp

deg:	node_deg_count.cpp $(HEADER_PATHS)/GDV_functions.hpp $(HEADER_PATHS)/bfs.hpp $(HEADER_PATHS)/print_disconnected_graph.hpp
	$(CXX) $(CXXFLAGS) -o node_deg $(ESSENS_PATHS) node_deg_count.cpp

stat:	static_assigner.cpp
	$(CXX) $(CXXFLAGS) -o static_assigner $(ESSENS_PATHS) static_assigner.cpp

cE:     create_edgelist.cpp
	$(CXX) -g -o cE.out $(ESSENS_PATHS) create_edgelist.cpp
sparse: sparsify.cpp
	$(CXX) -g  -o sparse.out $(ESSENS_PATHS) sparsify.cpp
tEx:    traversalEx.cpp
	$(CXX) -g  -o tEx.out $(ESSENS_PATHS) traversalEx.cpp
bfs:    BFS.cpp
	$(CXX) -g -o bfs.out $(ESSENS_PATHS) BFS.cpp
chChk:  chordalchk.cpp
	$(CXX) -g -o chChk.out $(ESSENS_PATHS) chordalchk.cpp
lpath:  longest_path.cpp
	$(CXX) -g -o lpath.out $(ESSENS_PATHS) longest_path.cpp
permF:  permute_file.cpp
	$(CXX) -g -o permF.out $(ESSENS_PATHS) permute_file.cpp
BFS_try: BFS.cpp
	$(CXX) -g -fopenmp -o BFST.out $(ESSENS_PATHS) BFS.cpp
test1:  test_SetOps.cpp
	$(CXX) $(BASIC_IO_INC) $(BASIC_SETOP_INC) $(BASIC_TRAVERSAL_INC) test_SetOps.cpp

clean:
	$(RM) *.out *.o graph_alignment node_deg static_assigner
