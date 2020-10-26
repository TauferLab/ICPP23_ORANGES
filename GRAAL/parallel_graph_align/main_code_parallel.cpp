// Required directives are to be added from ESSENS to use the functionality 

#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
#include "/home/pnbell/Src_GraphAlignment/GRAAL/headers/GDV_functions.hpp"
#include "/home/pnbell/Src_GraphAlignment/GRAAL/headers/class_definitions.hpp"
#include <time.h>
#include <mpi.h>

void Calculate_GDV(int ,A_Network ,vector<OrbitMetric>&, GDVMetric&);
void readin_orbits(  ifstream* ,vector<OrbitMetric>* );
void convert_string_vector_int(string* , vector<int>* ,string );
using namespace std;

int main(int argc, char *argv[]) {

  clock_t tStart = clock();
  clock_t q, q1, q2,t;
  GDV_functions gdvf;
  /* Accepts file as input. 
     Input should be in the format of ( node1 node2 weight ). 
     By default, weight should be 1
     Should be sorted. */ 
  ifstream the_file1 ( argv[1] ); 
  if (!the_file1.is_open() ) { 
    cout<<"INPUT ERROR:: Could not open the graph input file\n";
  }

   ifstream the_file2 ( argv[2] ); 
  if (!the_file2.is_open() ) { 
    cout<<"INPUT ERROR:: Could not open the orbit file\n";
  }
  
  vector<OrbitMetric> orbits;
  readin_orbits(&the_file2,&orbits);
  // Objects for testing orbit creation 
  // print_vector(orbits[1].orbitDegree);
  // vector<OrbitMetric> filter_o = gdvf.orbit_filter(&orbits,3);

  int numtasks, rank, dest, source, rc, count, tag=0;
  MPI_Status Stat;   // required variable for receive routines                                                                                                                                                                          

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Read in Network: Reads the file and converts it into a network of type A_Network*/
  A_Network X;
  readin_network(&X,argv[1],0,-1);  
  GDV_functions test_gdvf;

  // objects for testing orbit list
  // vector<OrbitMetric> filtered_orbits ;
  // gdvf.orbit_filter(orbits,3,filtered_orbits);
  // print_vector(filtered_orbits[1].orbitDistance);

  // Objects for testing GDV induced subgraph function
  //A_Network subgraph;
  //vector<int> subgraph_nodes;
  //subgraph_nodes.push_back(0);
  //subgraph_nodes.push_back(1);
  //subgraph_nodes.push_back(2);
  //subgraph_nodes.push_back(7);
  //gdvf.inducedSubgraph(X, subgraph_nodes, subgraph);
  //print_disconnected_network(subgraph);
  // cout<<"subgraph for 0,1,2,6"<<endl;
  // print_network(subgraph);

  // Objects for testing connectedness function
  //bool is_connected = false;
  //gdvf.isConnected(subgraph, is_connected);
  //if (is_connected) {
  //  cout << "subgraph connected" << endl;
  //} else {
  //  cout << "subgraph disconnected" << endl;
  //}

  // Objects for testing degree signature
  //vector<int> degree_sig;
  //test_gdvf.degree_signature(X, degree_sig);
  //print_vector(degree_sig);

  // Objects for testing distance signature
  //vector<int> distance_sig;
  //test_gdvf.distance_signature(2, X, distance_sig);
  //print_vector(distance_sig);
  
  // for (int i:X)
  // {
  //   // Calculate_GDV(i,X);
  // }

  //for (ADJ_Bundle node:X) {

  // Set up parallelization
  int comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  int nodes_per_proc = int(X.size())/comm_size;
  int assigned = nodes_per_proc * comm_size;
  int remain = X.size() - (nodes_per_proc * comm_size);
  
  // Define Number of nodes to send to each rank
  // 
  // int assigned_node_count = nodes_per_proc;
  // int remain_node_count = nodes_per_proc - 1;
  // define vector of needed count per process so each process as nodes_per_proc
  // add one to each of first remain entries of vector
  // have each process do number of nodes in vector entry of corresponding index

  // Assign count of nodes for each process
  vector<int> per_proc(comm_size, nodes_per_proc);
  int starting_node = 0;
  //per_proc.push_back(nodes_per_proc + 1);
  for (int i = 1; i < remain; i++) {
    per_proc[i] += 1;
  }
  for (int i = 0; i < rank; i++) {
    starting_node += per_proc[i];
  }

  // Have process run the number of nodes for itself
  int node = starting_node;
  for (int i = 0; i < per_proc[rank]; i++) { 
    vector<int> GDV_1;
    GDVMetric gdvMetric(node,GDV_1);
    Calculate_GDV(node,X,orbits,gdvMetric);
    cout<<"gdv for node "<<node<<endl;
    print_vector(gdvMetric.GDV);
    node += 1;
  }
  
  if (rank == 0) {
    // MPI_Gather from each relevant rank to make final calculations
  }

  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  
  //delete send_buff;
  //delete rec_buff;
  MPI_Finalize(); 
  return 0;

}

void Calculate_GDV(int node,A_Network Graph,vector<OrbitMetric> &orbits, GDVMetric &gdvMetric)
{
    GDV_functions gdvf;
    vector<int> gdv(orbits.size(),0);
    // printf("calculating GDV for node %d\n",node);
    vector<int> neighbours;
    gdvf.find_neighbours(node,Graph,4,&neighbours);
    //print_vector(neighbours);
    int set[neighbours.size()]; 
    std::copy( neighbours.begin(), neighbours.end(), set );
    int numElements = *(&set + 1) - set;
    for (int node_count = 1; node_count < 5; node_count++)
    {
      vector<vector<int>> combinationsList;
      gdvf.find_combinations(set, numElements,node_count,&combinationsList);
      // cout<<"Node count is "<<node_count<<endl;
      // cout<<"total combinations are : "<<combinationsList.size()<<endl;
      for (vector<int> combination : combinationsList)
      {
        A_Network induced_sgraph;
        vector<int> subgraph_degree_signature;
        vector<int> subgraph_distance_signature;
        bool is_connected = false;
        combination.push_back(node);
        gdvf.inducedSubgraph(Graph, combination, induced_sgraph);
        gdvf.isConnected(induced_sgraph, is_connected);
        if(is_connected)
        {
            gdvf.degree_signature(induced_sgraph,subgraph_degree_signature);
            gdvf.distance_signature(node,induced_sgraph,subgraph_distance_signature);
            vector<OrbitMetric> filter_orbits;
            gdvf.orbit_filter(orbits,node_count+1,filter_orbits);
            for(OrbitMetric orbit: filter_orbits)
            {
              sort(orbit.orbitDegree.begin(),orbit.orbitDegree.end());
              sort(subgraph_degree_signature.begin(),subgraph_degree_signature.end());
              if( orbit.orbitDistance == subgraph_distance_signature && 
                  orbit.orbitDegree == subgraph_degree_signature)
              {
                gdv[orbit.orbitNumber] +=1;
                break;
              }
            }
        }
      }
    }
    gdvMetric.GDV = gdv;
    gdvMetric.node = node;
}

//This method takes the file and converts it into orbits and saves in output
void readin_orbits(ifstream *file,vector<OrbitMetric>* output )
{
  string line;
  string signature_delimiter;
  string internal_delimiter;
  signature_delimiter = "/";
  internal_delimiter= ",";
  while(std::getline(*file,line))
  {
    string s= line;
    size_t pos = 0;
    vector<vector<int>> vector_line;
    do
    {
      vector<int> segment; 
      string token;
      pos = s.find(signature_delimiter);
      token = s.substr(0, pos);
      token.erase(remove(token.begin(), token.end(), '['), token.end());
      token.erase(remove(token.begin(), token.end(), ']'), token.end());
      convert_string_vector_int(&token,&segment,internal_delimiter);
      s.erase(0, pos + signature_delimiter.length());
      vector_line.push_back(segment);
    }
    while (pos!= std::string::npos);

    OrbitMetric orbMetric(vector_line[0][0],vector_line[1],vector_line[2]);
    output->push_back(orbMetric);
  }

}

//This method converts a string containing integers to a vector of integers
void convert_string_vector_int(string* str, vector<int>* output,string delimiter)
{
  size_t pos = 0;
  int token;
  string s;
  s = *str;
  do
  {
    pos = s.find(delimiter);
    token = stoi(s.substr(0, pos));
    output->push_back(token);
    s.erase(0, pos + delimiter.length());
  }
  while (pos != std::string::npos);

}



