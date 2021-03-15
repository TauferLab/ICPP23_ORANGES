// Required directives are to be added from ESSENS to use the functionality 

#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
#include "/home/pnbell/Src_GraphAlignment/GRAAL/headers/GDV_functions.hpp"
#include "/home/pnbell/Src_GraphAlignment/GRAAL/headers/class_definitions.hpp"
#include <time.h>
#include <stdlib.h>
#include <ctime>
#include <mpi.h>
#include <fstream>

#define MAX_COMM_SIZE 32
//#define DEBUG

void Calculate_GDV(int ,A_Network ,vector<OrbitMetric>&, GDVMetric&);
void readin_orbits(  ifstream* ,vector<OrbitMetric>* );
void convert_string_vector_int(string* , vector<int>* ,string );
void Similarity_Metric_calculation_for_two_graphs(A_Network, A_Network,vector<OrbitMetric>, string, string);
double GDV_distance_calculation(GDVMetric&, GDVMetric&);
void metric_formula(GDVMetric&, double*);
void GDV_vector_calculation(A_Network,vector<GDVMetric>*,  vector<OrbitMetric>, const char*);

// Define variables for keeping track of time for load imbalancing tests.
/*double total_time_taken_mpi[MAX_COMM_SIZE] = {};
double vec_calc_prior_gather_mpi[MAX_COMM_SIZE] = {};
double vec_calc_post_gather_mpi[MAX_COMM_SIZE] = {};
double gdv_calc_mpi[MAX_COMM_SIZE] = {};*/
double total_time_taken;
double vec_calc_prior_gather;
double vec_calc_post_gather;
int vec_calc_avg_node_deg;

using namespace std;

int main(int argc, char *argv[]) {

  //  if (rank == 0) {
  //clock_t out_tStart = clock();
  time_t now = time(0);
    //}
  clock_t tStart = clock();
  clock_t q, q1, q2,t;
  GDV_functions gdvf;
  /* Accepts file as input. 
     Input should be in the format of ( node1 node2 weight ). 
     By default, weight should be 1
     Should be sorted. */ 
  ifstream the_file0 ( argv[1] );
  if (!the_file0.is_open() ) {
    cout<<"INPUT ERROR:: Could not open the graph input file\n";
  }

  ifstream the_file1 ( argv[2] ); 
  if (!the_file1.is_open() ) { 
    cout<<"INPUT ERROR:: Could not open the graph input file\n";
  }

   ifstream the_file2 ( argv[3] ); 
  if (!the_file2.is_open() ) { 
    cout<<"INPUT ERROR:: Could not open the orbit file\n";
  }

  //ifstream the_file3 ( argv[4] + "time_results.txt" );
  //if (!the_file3.is_open() ) {
  //  cout << "INPUT ERROR:: Could not open the time recording file\n";
  //}

  vector<OrbitMetric> orbits;
  readin_orbits(&the_file2,&orbits);
  // Objects for testing orbit creation 
  // print_vector(orbits[1].orbitDegree);
  // vector<OrbitMetric> filter_o = gdvf.orbit_filter(&orbits,3);

  A_Network X;
  A_Network Y;
  readin_network(&X,argv[1],0,-1);
  readin_network(&Y,argv[2],0,-1);
  GDV_functions test_gdvf;

  //clock_t out_tStart = clock();

  int numtasks, rank, dest, source, rc, count, tag=0;
  MPI_Status Stat;   // required variable for receive routines                                                                                                                                                                          

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  cout << "MPI Setup" << endl;

  /* Read in Network: Reads the file and converts it into a network of type A_Network*/
  //A_Network X;
  //A_Network Y;
  //readin_network(&X,argv[1],0,-1);  
  //readin_network(&Y,argv[2],0,-1);
  //GDV_functions test_gdvf;
  //if (rank == 0) {
  //print_disconnected_network(X);
  //print_disconnected_network(Y);
  //}

  for (int i = 0; i < X.size(); i++) {
    for (int j = 0; j < X[i].ListW.size(); j++) {
      for (int k = 0; k < j; k++) {
	if (X[i].ListW[j].first == X[i].ListW[k].first) {
	  X[i].ListW.erase(X[i].ListW.begin() + j);
	}
      }
    }
  }
  for (int i = 0; i < Y.size(); i++) {
    for (int j = 0; j < Y[i].ListW.size(); j++) {
      for (int k = 0; k < j; k++) {
        if (Y[i].ListW[j].first == Y[i].ListW[k].first) {
          Y[i].ListW.erase(Y[i].ListW.begin() + j);
        }
      }
    }
  }

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


  // Get data on graph names
  string delimiter = "/";
  string graph_name1 = argv[1];
  string graph_name2 = argv[2];
  
  size_t pos1 = 0;
  string token1;
  while ((pos1 = graph_name1.find(delimiter)) != string::npos) {
    token1 = graph_name1.substr(0, pos1);
    //cout << token << endl;
    graph_name1.erase(0, pos1 + delimiter.length());
  }
  
  size_t pos2 = 0;
  string token2;
  while ((pos2 = graph_name2.find(delimiter)) != string::npos) {
    token2 = graph_name2.substr(0, pos2);
    graph_name2.erase(0, pos2 + delimiter.length());
  }

  #ifdef DEBUG
    cout << "Starting Similarity Metric Calculation on Rank: " << rank << endl;
  #endif

  // Perform Similarity Calculations
  Similarity_Metric_calculation_for_two_graphs(X,Y,orbits, graph_name1, graph_name2);

  #ifdef DEBUG
    cout << "Finished Similarity Metric Calculation on Rank: " << rank << endl;
  #endif

  // Set up for and Perform Runtime Management and Gather
  double* time_buff = NULL;
  int num_times = 3;
  double send_times[num_times];
  send_times[0] = total_time_taken;
  send_times[1] = vec_calc_prior_gather;
  send_times[2] = vec_calc_post_gather;
  if (rank == 0) {
    time_buff = (double *)malloc(numtasks*num_times*sizeof(double));
  }
  //cout << "Rank " << rank << " made it to last gather." << endl;
  MPI_Gather(send_times, num_times, MPI_DOUBLE, time_buff, num_times, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  #ifdef DEBUG
    cout << "Finished time gathering and prepped for time-keeping on Rank: " << rank << endl;
  #endif

  // Handle output of timing results
  if (rank == 0) {
    
    // Get date of run
    tm *ltm = localtime(&now);
    int year = 1900 + ltm->tm_year;
    int month = 1 + ltm->tm_mon;
    int day = ltm->tm_mday;

    // File IO to Record Run Data
    // Date \t n_procs \t graph_name1 \t graph_name2 \t n_nodes \t runtime(s) 
    ofstream myfile;
    myfile.open(argv[4], ios_base::app);
    if (!myfile.is_open() ) {
      cout << "INPUT ERROR:: Could not open the time recording file\n";   
    }
  
    if (myfile.is_open()) {
      myfile << month << "/" << day << "/" << year << "\t" << numtasks << "\t" << graph_name1 << "\t" << graph_name2 << "\t" << X.size() << "\t\t" << total_time_taken << " \n";
      myfile.close();
    }
    else { 
      cout << "Out File Did Not Open";
    }

    // Print out rank specific runtime data
    string time_file = "runtimes_rec.txt";
    myfile.open(time_file, ofstream::trunc);
    if (!myfile.is_open()) {
      cout << "INPUT ERROR:: Could not open the local time recording file\n";
    }
    myfile << "Time Taken in Similarity Metric Calculation = " << " \n";
    for (int i = 0; i < numtasks; i++) {
      myfile << i << " " << time_buff[num_times*i+1] << " " << time_buff[num_times*i+2] << " " << time_buff[num_times*i] << " \n";
    }
    myfile.close();

  }


  //printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  cout << "Time taken on rank " << rank << " = " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
  //MPI_Barrier( MPI_COMM_WORLD );
  MPI_Finalize(); 
  return 0;

}

void Similarity_Metric_calculation_for_two_graphs(A_Network graph1, A_Network graph2,vector<OrbitMetric> orbits, string graph_tag1, string graph_tag2)
{

  //clock_t out_tStart = clock();
  MPI_Barrier( MPI_COMM_WORLD );
  double out_tStart = MPI_Wtime();

  vector<GDVMetric> graph1_GDV;
  vector<GDVMetric> graph2_GDV;
  int rankm, numtasksm;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankm);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasksm);

  GDV_vector_calculation(graph1, &graph1_GDV, orbits, "graph1"); 
  //cout << "Rank " << rankm << " finished graph1" << endl;
  GDV_vector_calculation(graph2, &graph2_GDV, orbits, "graph2"); 
  //cout << "Finished Vector Calculations" << endl;

  // Calculate Similarities in GDVs
  int m = (int)graph1_GDV.size();
  int n = (int)graph2_GDV.size();
  double sim_mat[m][n];
  for (GDVMetric &gdvm1: graph1_GDV) {
    for(GDVMetric &gdvm2: graph2_GDV) {
      sim_mat[gdvm1.node][gdvm2.node]= GDV_distance_calculation(gdvm1,gdvm2);
    }
  }


  // Measure final total time 
  total_time_taken = (double)(MPI_Wtime() - out_tStart);

  #ifdef DEBUG
    cout << "Finished similarity matrix calculation and prepped for fileIO on Rank: " << rankm << endl;
  #endif

  // Perform File I/O for Output Data
  if (rankm == 0) {

    // Start by Printing out Similarity Matrix into a file
    ofstream myfile; 
    string filename = "out_similarity_matrix.txt";
    myfile.open(filename);
    for(int i=1; i<m;i++) {
      for(int j=1;j<n;j++) {
	myfile<<" { "<<sim_mat[i][j]<<" } ";
      }
      myfile<<"||"<<endl;
    }
    myfile.close();
    
    
    // Print out GDVs into files
    string gdv_file1 = "out_gdv_1_" + graph_tag1 + ".txt";
    string gdv_file2 = "out_gdv_2_" + graph_tag2 + ".txt";
    myfile.open(gdv_file1, ofstream::trunc);
    for (int i = 0; i < graph1_GDV.size(); i++) {
      for (int j = 0; j< graph1_GDV[i].GDV.size(); j++) {
	myfile << graph1_GDV[i].GDV[j] << " ";
      }
      myfile << endl;
    }
    myfile.close();
    myfile.open(gdv_file2, ofstream::trunc);
    for (int i = 0; i < graph2_GDV.size(); i++) {
      for (int j = 0; j< graph2_GDV[i].GDV.size(); j++) {
        myfile << graph2_GDV[i].GDV[j] << " ";
      }
      myfile << endl;
    }
    myfile.close();

  }
}

double GDV_distance_calculation(GDVMetric &gdvm1, GDVMetric &gdvm2)
{   
  double gdv1_score;
  double gdv2_score;
  metric_formula(gdvm1,&gdv1_score);
  metric_formula(gdvm2,&gdv2_score);
  
  double similarity_score = abs(gdv1_score - gdv2_score);

  return similarity_score;
}

void metric_formula(GDVMetric &gdvm, double* gdv_score)
{
  int sum=0;
  //  Formula used here is the vector norm 2 or the l2 norm. 
  //  We square each value and do summation. Now take a square root.
  for(int x: gdvm.GDV)
    {
      int number;
      number = x*x;
      sum = sum + number;
    }
  *gdv_score = sqrt(sum);
}

void GDV_vector_calculation(A_Network graph,vector<GDVMetric>* graph_GDV,  vector<OrbitMetric> orbits, const char* graph_name)
{

  double vec_calc_start = MPI_Wtime();

  //GDV_functions gdvfunc;
  //bool connected = false;
  //gdvfunc.isConnected(graph, connected);
  //cout << "graph connectedness: " << connected << endl;

  // Set up parallelization               
  int comm_size, rankn;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankn);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  int nodes_per_proc = int(graph.size())/comm_size;
  int assigned = nodes_per_proc * comm_size;
  int remain = graph.size() - (nodes_per_proc * comm_size);
  graph_GDV->clear();


  // Assign count of nodes for each process                        
  vector<int> per_proc(comm_size, nodes_per_proc);
  int starting_node = 0;
  for (int i = 0; i < remain; i++) {
    per_proc[i] += 1;
  }
  for (int i = 0; i < rankn; i++) {
    starting_node += per_proc[i];
  }

  #ifdef DEBUG
    cout << "Assigned Nodes for GDV Calculation on Rank: " << rankn << endl;
  #endif
 
  // Have process run the number of nodes for itself     
  int node = starting_node;
  vector<GDVMetric> temp_gdvs;
  int gdv_length;
  int per = per_proc[rankn];
  for (int i = 0; i < per; i++) {
    vector<int> GDV_1;
    GDVMetric gdvMetric(node,GDV_1);
    Calculate_GDV(node,graph,orbits,gdvMetric);
    if (i == 0) {gdv_length = gdvMetric.GDV.size();}
    temp_gdvs.push_back(gdvMetric);
    node += 1;
  }

  #ifdef DEBUG
    cout << "Finished GDV Calculation on Rank: " << rankn << endl;
  #endif

  // Set up displacement array for gatherv
  int per_displ[comm_size];
  per_displ[0] = 0;
  int displ_sum = 0;
  for (int i = 1; i < comm_size; i++) {
    displ_sum += per_proc[i-1];
    per_displ[i] = displ_sum;
  }

  // Set up receive counts for each rank in an array for gatherv
  int rcv_sizes[comm_size];
  int graph_size = graph.size();
  for (int i = 0; i < comm_size; i++) {
    rcv_sizes[i] = per_proc[i];
  }

  // Set up the nodes and gdv data to send using an array for gatherv
  int send_gdvs[per * gdv_length];
  int send_nodes[per];
  for (int i = 0; i < temp_gdvs.size(); i++) {
    send_nodes[i] = temp_gdvs[i].node;
    for (int j = 0; j < gdv_length; j++) {
      send_gdvs[i*gdv_length+j] = temp_gdvs[i].GDV[j];
    }
  }


  // Declare and allocate space for receive buffers for gather
  int gdv_size = sizeof(GDVMetric);
  int *gdvs = NULL;
  int *nodes = NULL;
  if (rankn == 0) {
    gdvs = (int *)calloc(gdv_length * graph.size(), sizeof(int));
    nodes = (int *)calloc(graph.size(), sizeof(int));
  }

  #ifdef DEBUG
    cout << "Prepped for MPI Gatherv on Rank: " << rankn << endl;
  #endif

  // Prepare for and implement gatherv to get gdvs from each rank
  vec_calc_prior_gather = MPI_Wtime() - vec_calc_start + vec_calc_prior_gather;
  //if (rankn == 0) {
  //  cout << "Rank 0 made it to gdv gathers." << endl;
  //}
  MPI_Gatherv(send_nodes, per, MPI_INT, nodes, rcv_sizes, per_displ, MPI_INT, 0, MPI_COMM_WORLD);
  for (int i = 0; i < comm_size; i++) {
    per_displ[i] = per_displ[i] * gdv_length;
    rcv_sizes[i] = rcv_sizes[i] * gdv_length;
  }
  MPI_Gatherv(send_gdvs, gdv_length * per, MPI_INT, gdvs, rcv_sizes, per_displ, MPI_INT, 0, MPI_COMM_WORLD);

  #ifdef DEBUG
    cout << "Finished MPI Gatherv on Rank: " << rankn << endl;
  #endif

  // Post processing of MPI gather
  if (rankn == 0) {
    
    // Get data from receive buffers into output GDV list
    vector<int> gdv();
    for (int i = 0; i < graph.size(); i++) {
      vector<int> gdv(&gdvs[i * gdv_length], &gdvs[i * gdv_length] + gdv_length);
      GDVMetric gdv_extract(nodes[i], gdv);
      graph_GDV->push_back(gdv_extract);
    }
   
    // Free dynamically allocated recieve buffers
    if (gdvs == NULL) {
      cout << "gdvs is null" << endl;
    } else {
    free(gdvs);
    } 
    if (nodes == NULL) {
      cout << "nodes is null" << endl;
    } else {
    free(nodes);
    }

    //cout << "Rank 0 made it to gdv gathers." << endl;
  }

  vec_calc_post_gather = MPI_Wtime() - vec_calc_start + vec_calc_post_gather;

  #ifdef DEBUG
    cout << "Finished GDV Vector Calc on Rank: " << rankn << endl;
  #endif

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
      for (vector<int> &combination : combinationsList)
      {
        A_Network induced_sgraph;
        vector<int> subgraph_degree_signature;
        vector<int> subgraph_distance_signature;
        bool is_connected = false;
        combination.push_back(node);
        gdvf.inducedSubgraph(Graph, combination, induced_sgraph);
	//cout << "Made subgraph" << endl;
        gdvf.isConnected(induced_sgraph, is_connected);
	//cout << "Finished connectedness check" << endl;
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



