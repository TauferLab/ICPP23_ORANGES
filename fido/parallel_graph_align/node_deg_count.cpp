#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
#include <stdlib.h>
#include <ctime>
#include <mpi.h>
#include <fstream>


float avg_deg;
int deg_count;


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


  MPI_Init(&argc,&argv);
  int comm_size, rankn;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankn);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);


/*  int nodes_per_proc = int(X.size())/comm_size;
  int assigned = nodes_per_proc * comm_size;
  int remain = X.size() - (nodes_per_proc * comm_size);
  vector<int> per_proc(comm_size, nodes_per_proc);
  int starting_node = 0;
  for (int i = 0; i < remain; i++) {
    per_proc[i] += 1;
  }
  for (int i = 0; i < rankn; i++) {
    starting_node += per_proc[i];
  } */
  int part_assignment;
  int node_counter = 0;
  vector<int> per_proc(comm_size, 0);
  vector<int> assigned_nodes;
  ifstream partfile;
  partfile.open(argv[3]);
  if (!partfile.is_open() ) {
    cout << "INPUT ERROR:: Could not open the partition file\n";
  }
  if (partfile.is_open()) {
    while (true) {
      partfile >> part_assignment;
      //#ifdef PARTITION
      //  cout << "Node " << node_counter << " has been assigned to rank " << part_assignment << endl;
      //#endif
      if (partfile.eof()) {
        break;
      }
      if (part_assignment == rankn) {
        assigned_nodes.push_back(X[node_counter].Row);
      }
      per_proc[part_assignment] += 1;
      node_counter += 1;
    }
    partfile.close();
  }
  

  // Calculate the average degree of nodes assigned from graph 1
  //int node = starting_node;
  int gdv_length;
  int temp_deg;
  int per = per_proc[rankn];
  for (int i = 0; i < assigned_nodes.size(); i++) {
    temp_deg = 0;
    deg_count += 1;
    temp_deg += X[assigned_nodes[i]].ListW.size();
    for (int j = 0; j < X.size(); j++) {
      for (int k = 0; k < X[j].ListW.size(); k++) {
	if (X[j].ListW[k].first == assigned_nodes[i]) {
	  temp_deg += 1;
	}
      }
    }
    avg_deg += temp_deg;
  }


/*  nodes_per_proc = int(Y.size())/comm_size;
  assigned = nodes_per_proc * comm_size;
  remain = Y.size() - (nodes_per_proc * comm_size);
  vector<int> per_proc2(comm_size, nodes_per_proc);
  starting_node = 0;
  for (int i = 0; i < remain; i++) {
    per_proc2[i] += 1;
  }
  for (int i = 0; i < rankn; i++) {
    starting_node += per_proc2[i];
  } */
  //int part_assignment;
  node_counter = 0;
  vector<int> per_proc2(comm_size, 0);
  vector<int> assigned_nodes2;
  //ifstream partfile;
  partfile.open(argv[4]);
  if (!partfile.is_open() ) {
    cout << "INPUT ERROR:: Could not open the partition file\n";
  }
  if (partfile.is_open()) {
    while (true) {
      partfile >> part_assignment;
      //#ifdef PARTITION
      //  cout << "Node " << node_counter << " has been assigned to rank " << part_assignment << endl;
      //#endif
      if (partfile.eof()) {
        break;
      }
      if (part_assignment == rankn) {
        assigned_nodes2.push_back(Y[node_counter].Row);
      }
      per_proc2[part_assignment] += 1;
      node_counter += 1;
    }
    partfile.close();
  }


  //node = starting_node;
  //temp_deg;
  per = per_proc2[rankn];
  for (int i = 0; i < assigned_nodes2.size(); i++) {
    temp_deg = 0;
    deg_count += 1;
    temp_deg += Y[assigned_nodes2[i]].ListW.size();
    for (int j = 0; j < Y.size(); j++) {
      for (int k = 0; k < Y[j].ListW.size(); k++) {
        if (Y[j].ListW[k].first == assigned_nodes2[i]) {
          temp_deg += 1;
        }
      }
    }
    avg_deg += temp_deg;
    //node += 1;
  }


  // Calculate Average Degree
  float temp_avg = avg_deg / deg_count;
  avg_deg = temp_avg;


  ofstream fileio;
  string filename = "deg_count_per_proc.txt";
  fileio.open(filename, ios_base::app);
  //cout << rankn << " " << avg_deg << " \n";
  fileio << rankn << " " << avg_deg << " " << deg_count << " \n";
  fileio.close();


  MPI_Finalize();

  return 0;

}










