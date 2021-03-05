#Change thos line to a different location if needed
#ESSENS =/mnt/c/users/lohit/desktop/RA/Src_GraphAlignment/ESSENS-master
ESSENS=/home/pnbell/Src_GraphAlignment/ESSENS-master
#ESSENS= ./ESSENS

#Use this to run on Linux
CXX=g++
CXXFLAGS=-g -fopenmp -O3

#Changed CXX flags to run on MAC
#CXX = /usr/local/opt/llvm/bin/clang++
#CXXFLAGS = -I/usr/local/opt/llvm/include -fopenmp
# LDFLAGS = -L/usr/local/opt/llvm/lib

# CXX = mpicxx
#CXXFLAGS = -openmp

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

all:    main_code.cpp headers/GDV_functions.hpp headers/bfs.hpp headers/print_disconnected_graph.hpp
	echo $(ESSENS_PATHS)
	$(CXX) $(CXXFLAGS) -o graal $(ESSENS_PATHS) main_code.cpp

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
	$(RM) *.out *.o graal
