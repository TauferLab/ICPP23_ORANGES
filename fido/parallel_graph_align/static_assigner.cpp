#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
#include <stdlib.h>
#include <ctime>
#include <mpi.h>
#include <fstream>


int main(int argc, char *argv[]) {


  cout << "Entered Main" << endl;
  ifstream the_file0 ( argv[1] );
  if (!the_file0.is_open() ) {
    cout<<"INPUT ERROR:: Could not open the graph input file\n";
  }

  ifstream the_file1 ( argv[2] );
  if (!the_file1.is_open() ) {
    cout<<"INPUT ERROR:: Could not open the graph input file\n";
  }


  A_Network X;
  A_Network Y;
  readin_network(&X,argv[1],0,-1);
  readin_network(&Y,argv[2],0,-1);
  cout << "Read Files" << endl;

  int assigned_nodes_X[X.size()] = {0};
  int assigned_nodes_Y[Y.size()] = {0};
  string assigner_style = argv[3];

  MPI_Init(&argc,&argv);
  int comm_size, rankn;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankn);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);


  if ( assigner_style == "mod_divide" ) {

    int node_counter = 0;
    int nodes_per_proc = int(X.size())/comm_size;
    int remain = X.size() - (nodes_per_proc * comm_size);
    vector<int> per_proc_x(comm_size, nodes_per_proc);
    for (int i = 0; i < remain; i++) {
      per_proc_x[i] += 1;
    }  

    //int assigned_nodes_X[X.size()] = {0};
    for (int i = 0; i < per_proc_x.size(); i++) {
      for (int j = 0; j < per_proc_x[i]; j++) {
	assigned_nodes_X[node_counter] = i;
	node_counter += 1;
      }
    }

    node_counter = 0;
    nodes_per_proc = int(X.size())/comm_size;
    remain = X.size() - (nodes_per_proc * comm_size);
    vector<int> per_proc_y(comm_size, nodes_per_proc);
    for (int i = 0; i < remain; i++) {
      per_proc_y[i] += 1;
    }                                
  
    //int assigned_nodes_Y[Y.size()] = {0};
    for (int i = 0; i < per_proc_y.size(); i++) {
      for (int j = 0; j < per_proc_y[i]; j++) {
        assigned_nodes_Y[node_counter] = i;
        node_counter += 1;
      }
    }    

  }

  if ( assigner_style == "round_robin" ) {
    
  }


  ofstream myfile;
  myfile.open(argv[4], ios_base::trunc);
  if (!myfile.is_open() ) {
    cout << "INPUT ERROR:: Could not open the partition file for graph 1\n";
  }

  if (myfile.is_open()) {
    for (int i = 0; i < X.size(); i++) {
      myfile << assigned_nodes_X[i] << "\n"; 
    }
    myfile.close();
  }

  myfile.open(argv[5], ios_base::trunc);
  if (!myfile.is_open() ) {
    cout << "INPUT ERROR:: Could not open the partition file for graph 2\n";
  }

  if (myfile.is_open()) {
    for (int i = 0; i < Y.size(); i++) {
      myfile << assigned_nodes_Y[i] << "\n";
    }
    myfile.close();
  }

  MPI_Finalize();

  return 0;

}








