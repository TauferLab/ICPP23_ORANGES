// Required directives are to be added from ESSENS to use the functionality 

#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
#include "GDV_functions.hpp"
#include "class_definitions.hpp"
#include "combinations.hpp"
#include "kokkos_functions.hpp"
#include "fido_compile_config.h"
#include <time.h>
#include <stdlib.h>
#include <ctime>
#include <mpi.h>
#include <fstream>

#define MAX_COMM_SIZE 32
#define GDV_LENGTH 22
#define CHUNK_SIZE 1
#define SIM_MAT_CUTOFF 1000
//#define DEBUG
//if (TRACK_RUNTIME) {
#define RUNTIME TRACK_RUNTIME
//}
#define NUM_THREADS 1


void Calculate_GDV(int ,A_Network ,vector<OrbitMetric>&, GDVMetric&);
void readin_orbits(  ifstream* ,vector<OrbitMetric>* );
void convert_string_vector_int(string* , vector<int>* ,string );
void Similarity_Metric_calculation_for_two_graphs(A_Network, A_Network,vector<OrbitMetric>, string, string);
double GDV_distance_calculation(GDVMetric&, GDVMetric&);
void metric_formula(GDVMetric&, double*);
void GDV_vector_calculation(A_Network,vector<GDVMetric>*,  vector<OrbitMetric>, const char*, int);

struct node_id_order {
  inline bool operator() (const GDVMetric& gdv_met_1, const GDVMetric& gdv_met_2) {
    return (gdv_met_1.node < gdv_met_2.node);
  }
};

// Define variables for keeping track of time for load imbalancing tests.
#if RUNTIME
double total_time_taken;
double vec_calc_communication_time[2] = {};
double pre_process_time;
double report_output_time;
double* vec_calc_computation_time_X;
double* vec_calc_computation_time_Y;
int* vec_calc_proc_assign_X;
int* vec_calc_proc_assign_Y;
#endif
//double vec_calc_prior_gather;
//double vec_calc_post_gather;
//int vec_calc_avg_node_deg;

using namespace std;

int main(int argc, char *argv[]) {


//  CombinationGenerator generator(5,3);
//  do {
//    printf("Combination: ");
//    auto& combo = generator.get_combination();
//    for(int i=0; i<generator.k; i++) {
//      printf("%d ", combo[i]);
//    }
//    printf("\n");
//  } while(!generator.next());
  MPI_Init(&argc,&argv);
//  Kokkos::initialize(argc, argv);
//  {

  //  if (rank == 0) {
  //clock_t out_tStart = clock();
  time_t now = time(0);
  #if RUNTIME
  double pre_tStart = MPI_Wtime();
  #endif
    //}
  //clock_t tStart = clock();
  //clock_t q, q1, q2,t;
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

//  ifstream the_file3 ( argv[4] + "time_results.txt" );
//  if (!the_file3.is_open() ) {
//    cout << "INPUT ERROR:: Could not open the time recording file\n";
//  }

  vector<OrbitMetric> orbits;
  readin_orbits(&the_file2,&orbits);
  // Objects for testing orbit creation 
  // print_vector(orbits[1].orbitDegree);
  // vector<OrbitMetric> filter_o = gdvf.orbit_filter(&orbits,3);

  //A_Network X;
  //A_Network Y;
  //readin_network(&X,argv[1],0,-1);
  //readin_network(&Y,argv[2],0,-1);
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

  // Allocate space for time recording by graph node

#if RUNTIME
  vec_calc_computation_time_X = new double[graphx.numRows()]();
  vec_calc_computation_time_Y = new double[graphy.numRows()]();
  vec_calc_proc_assign_X = new int[graphx.numRows()]();
  vec_calc_proc_assign_Y = new int[graphy.numRows()]();
#endif

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

  //vec_calc_communication_time = 0;
  //vec_calc_computation_time = 0;

  #if RUNTIME
  pre_process_time = (double)(MPI_Wtime() - pre_tStart);
  #endif
  // Perform Similarity Calculations
  Similarity_Metric_calculation_for_two_graphs(X,Y,orbits, graph_name1, graph_name2);

  #ifdef DEBUG
    cout << "Finished Similarity Metric Calculation on Rank: " << rank << endl;
  #endif

  // Set up for and Perform Runtime Management and Gather
  double* time_buff = NULL;
  int num_times = 5;
  double send_times[num_times];
//  if (rank == 0) {
//    vec_calc_computation_time = 0;
//  }
//  cout << "Runtime on rank " << rank << " for graph 1 = " << vec_calc_communication_time[0] << endl;
//  cout << "Runtime on rank " << rank << " for graph 2 = " << vec_calc_communication_time[1] << endl;
#if RUNTIME
  send_times[0] = total_time_taken;
  send_times[1] = vec_calc_communication_time[0];
  send_times[2] = vec_calc_communication_time[1];
  send_times[3] = pre_process_time;
  send_times[4] = report_output_time;
  if (rank == 0) {
    time_buff = (double *)malloc(numtasks*num_times*sizeof(double));
  }
  //cout << "Rank " << rank << " made it to last gather." << endl;
  MPI_Gather(send_times, num_times, MPI_DOUBLE, time_buff, num_times, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  #ifdef DEBUG
    cout << "Finished time gathering and prepped for time-keeping on Rank: " << rank << endl;
  #endif

  // Handle output of timing results
  #if RUNTIME
  if (rank == 0) {
    
    // Get date of run
    tm *ltm = localtime(&now);
    int year = 1900 + ltm->tm_year;
    int month = 1 + ltm->tm_mon;
    int day = ltm->tm_mday;

    // File IO to Record Run Data
    // Date \t n_procs \t graph_name1 \t graph_name2 \t n_nodes \t runtime(s) 
    ofstream myfile;
    /*myfile.open(argv[4], ios_base::app);
    if (!myfile.is_open() ) {
      cout << "INPUT ERROR:: Could not open the time recording file\n";   
    }
  
    if (myfile.is_open()) {
      myfile << month << "/" << day << "/" << year << "\t" << numtasks << "\t" << graph_name1 << "\t" << graph_name2 << "\t" << X.size() << "\t\t" << total_time_taken << " \n";
      myfile.close();
    }
    else { 
      cout << "Out File Did Not Open";
    }*/

    // Print out rank specific runtime data
    string time_file = "runtime_data/runtimes_rec_over.txt";
    myfile.open(time_file, ofstream::trunc);
    if (!myfile.is_open()) {
      cout << "INPUT ERROR:: Could not open the local time recording file\n";
    }
    myfile << "Time Taken in Similarity Metric Calculation = " << " \n";

    myfile << "Rank\tPre-Process\tReport Data\tGraph 1\tGraph 2\tTotal\n";
    for (int i = 0; i < numtasks; i++) {
      myfile << i << " " << time_buff[num_times*i+3] << " " << time_buff[num_times*i+4] << " " << time_buff[num_times*i+1] << " " << time_buff[num_times*i+2] << " " << time_buff[num_times*i] << " \n";
    }
    myfile.close();

  }

  string computation_time_file_x = "runtime_data/runtimes_rec_x_" + to_string(rank) + ".txt";
  string computation_time_file_y = "runtime_data/runtimes_rec_y_" + to_string(rank) + ".txt";
  ofstream myfile;
  myfile.open(computation_time_file_x, ofstream::trunc);
  if (!myfile.is_open()) {
    cout << "INPUT ERROR:: Could not open the local time recording file\n";
  }
  if (myfile.is_open()) {
    for (int i = 0; i < X.size(); i++) {
      myfile << X[i].Row << " " << vec_calc_proc_assign_X[i] << " " << vec_calc_computation_time_X[i] << endl;
    }
    myfile.close();
  }
  myfile.open(computation_time_file_y, ofstream::trunc);
  if (!myfile.is_open()) {
    cout << "INPUT ERROR:: Could not open the local time recording file\n";
  }
  if (myfile.is_open()) {
    for (int i = 0; i < Y.size(); i++) {
      myfile << Y[i].Row << " " << vec_calc_proc_assign_Y[i] << " " << vec_calc_computation_time_Y[i] << endl;
    }
    myfile.close();
  }
  #endif

  //printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  //cout << "Time taken on rank " << rank << " = " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
  //MPI_Barrier( MPI_COMM_WORLD );
  #if RUNTIME
  delete[] vec_calc_computation_time_X;
  delete[] vec_calc_computation_time_Y;
  delete[] vec_calc_proc_assign_X;
  delete[] vec_calc_proc_assign_Y;

  MPI_Finalize(); 
  return 0;

}

void Similarity_Metric_calculation_for_two_graphs(A_Network graph1, A_Network graph2,vector<OrbitMetric> orbits, string graph_tag1, string graph_tag2)
{

  //clock_t out_tStart = clock();
  MPI_Barrier( MPI_COMM_WORLD );
  #ifdef RUNTIME
  double out_tStart = MPI_Wtime();
  #endif

  vector<GDVMetric> graph1_GDV;
  vector<GDVMetric> graph2_GDV;
  int rankm, numtasksm;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankm);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasksm);

  int graph_counter = 1;
  GDV_vector_calculation(graph1, &graph1_GDV, orbits, "graph1", graph_counter); 

//if(rankm == 0) {
//printf("GDV node 1\n");
//for(int i=0; i<graph1_GDV[1].GDV.size(); i++) {
//  printf("%d, ", graph1_GDV[1].GDV[i]);
//}
//printf("\n");
//printf("GDV node 10\n");
//for(int i=0; i<graph1_GDV[10].GDV.size(); i++) {
//  printf("%d, ", graph1_GDV[10].GDV[i]);
//}
//printf("\n");
//}
  //cout << "Rank " << rankm << " finished graph1" << endl;
//printf("Rank %d Finished graph 1\n", rankm);
  graph_counter = 2;
  GDV_vector_calculation(graph2, &graph2_GDV, orbits, "graph2", graph_counter); 
//printf("Finished graph 2\n");
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
  #ifdef RUTNIME
  total_time_taken = (double)(MPI_Wtime() - out_tStart);
  #endif

  #ifdef DEBUG
    cout << "Finished similarity matrix calculation and prepped for fileIO on Rank: " << rankm << endl;
  #endif

  // Perform File I/O for Output Data
  if (rankm == 0) {

    // Start by Printing out Similarity Matrix into a file
    ofstream myfile; 
    if ( (m > SIM_MAT_CUTOFF) || (n > SIM_MAT_CUTOFF) ) {
    string filename = "out_similarity_matrix.txt";
    myfile.open(filename);
    for(int i=1; i<m;i++) {
      for(int j=1;j<n;j++) {
	myfile<<" { "<<sim_mat[i][j]<<" } ";
      }
      myfile<<"||"<<endl;
    }
    myfile.close();
    }
    
    
    // Print out GDVs into files
    string gdv_file1 = "out_gdv_1_" + graph_tag1;
    string gdv_file2 = "out_gdv_2_" + graph_tag2;
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

/*KOKKOS_INLINE_FUNCTION void 
kokkos_Similarity_Metric_calculation_for_two_graphs(const matrix_type& graph1, 
                                                    const matrix_type& graph2, 
                                                    const Orbits& orbits, 
                                                    string graph_tag1, 
                                                    string graph_tag2) {
  //clock_t out_tStart = clock();
  #if RUNTIME
  MPI_Barrier( MPI_COMM_WORLD );
  double out_tStart = MPI_Wtime();
  #endif

  int num_orbits = orbits.degree.extent(0);
//cout << "Number of orbits: " << num_orbits << ", number of vertices: " << graph1.numRows() << endl;
  GDVs graph1_GDV("GDV for graph1", graph1.numRows(), num_orbits);
  GDVs graph2_GDV("GDV for graph2", graph2.numRows(), num_orbits);

  int rankm, numtasksm;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankm);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasksm);

  int graph_counter = 1;
  kokkos_GDV_vector_calculation(graph1, graph1_GDV, orbits, "graph1", graph_counter); 
  MPI_Barrier(MPI_COMM_WORLD);

if(rankm == 0)
  cout << "Done with GDV calculation for graph 1\n";

//if(rankm == 0) {
//cout << "Post kokkos_GDV_vector_calculation: graph 1" << endl;
//for(int i=0; i<graph1_GDV.extent(0); i++) {
//  cout << "Node " << i << ": ";
//  for(int j=0; j<graph1_GDV.extent(1); j++) {
//    cout << graph1_GDV(i,j) << " ";
//  }
//  cout << endl;
//}
//cout << endl;
//}

  graph_counter = 2;
  kokkos_GDV_vector_calculation(graph2, graph2_GDV, orbits, "graph2", graph_counter); 
  MPI_Barrier(MPI_COMM_WORLD);

if(rankm == 0)
  cout << "Done with GDV calculation for graph 2\n";

//if(rankm == 0) {
//cout << "Post kokkos_GDV_vector_calculation: graph 2" << endl;
//for(int i=0; i<graph2_GDV.extent(0); i++) {
//  cout << "Node " << i << ": ";
//  for(int j=0; j<graph2_GDV.extent(1); j++) {
//    cout << graph2_GDV(i,j) << " ";
//  }
//  cout << endl;
//}
//cout << endl;
//}

  int m = graph1_GDV.extent(0);
  int n = graph2_GDV.extent(0);
  Kokkos::View<double**> sim_mat("Similarity matrix", m, n);
if(rankm == 0) {
  Kokkos::MDRangePolicy<Kokkos::Rank<2>> mdrange({0,0}, {m,n});
  Kokkos::parallel_for("Compute sim_mat", mdrange, KOKKOS_LAMBDA(const int i, const int j) {
      auto g1_subview = Kokkos::subview(graph1_GDV, i, Kokkos::ALL);
      auto g2_subview = Kokkos::subview(graph2_GDV, j, Kokkos::ALL);
      sim_mat(i,j) = kokkos_GDV_distance_calculation(g1_subview, g2_subview);
  });
}
if(rankm == 0)
  cout << "Created similarity matrix\n";
//cout << "Similarity matrix\n";
//for(int i=0; i<m; i++) {
//  for(int j=0; j<n; j++) {
//    cout << sim_mat(i,j) << " ";
//  }
//  cout << endl;
//}
//cout << endl;
//  });

  // Measure final total time 
  #if RUNTIME
  total_time_taken = (double)(MPI_Wtime() - out_tStart);
  #endif

  #ifdef DEBUG
    cout << "Finished similarity matrix calculation and prepped for fileIO on Rank: " << rankm << endl;
  #endif

  double report_tStart = MPI_Wtime();

  // Perform File I/O for Output Data
  if (rankm == 0) {

    // Start by Printing out Similarity Matrix into a file
    ofstream myfile; 
    if ( (m > SIM_MAT_CUTOFF) || (n > SIM_MAT_CUTOFF) ) {
    string filename = "out_similarity_matrix.txt";
    myfile.open(filename);
    for(int i=1; i<m;i++) {
      for(int j=1;j<n;j++) {
        myfile<<" { "<<sim_mat(i,j)<<" } ";
      }
      myfile<<"||"<<endl;
    }
    myfile.close();
    }
    
    
    // Print out GDVs into files
    string gdv_file1 = "out_gdv_1_" + graph_tag1 + ".txt";
    string gdv_file2 = "out_gdv_2_" + graph_tag2 + ".txt";
    myfile.open(gdv_file1, ofstream::trunc);
    for (int i = 0; i < graph1_GDV.extent(0); i++) {
      for (int j = 0; j< graph1_GDV.extent(1); j++) {
        myfile << graph1_GDV(i,j) << " ";
      }
      myfile << endl;
    }
    myfile.close();
    myfile.open(gdv_file2, ofstream::trunc);
    for (int i = 0; i < graph2_GDV.extent(0); i++) {
      for (int j = 0; j< graph2_GDV.extent(1); j++) {
        myfile << graph2_GDV(i,j) << " ";
      }
      myfile << endl;
    }
    myfile.close();

  }
}*/
  #if RUNTIME
  report_output_time = (double)(MPI_Wtime() - report_tStart);
  #endif
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

void GDV_vector_calculation(A_Network graph,vector<GDVMetric>* graph_GDV,  vector<OrbitMetric> orbits, const char* graph_name, int graph_counter)
{

  // Set up parallelization               
  int comm_size, rankn;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankn);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  int tag = 11;
  int graph_size = graph.size();
  graph_GDV->clear();

//  int max_node = 0;
//  for(int i=0; i<graph.size(); i++) {
//    if(max_node < graph[i].Row)
//      max_node = graph[i].Row;
//  }
//  if(max_node < graph.size()) {
//    graph_GDV->resize(max_node+1);
//  } else {
//    graph_GDV->resize(max_node+2);
//  }
  for(int i=0; i<graph.size(); i++) {
    vector<int> gdv_vec(GDV_LENGTH, 0);
    graph_GDV->push_back(GDVMetric(graph[i].Row, gdv_vec));
  }
  //graph_counter += 1;
  //cout << "Graph counter on rank " << rankn << " = " << graph_counter << endl;

#ifdef RUNTIME
  double process_ends_communication;
  double vec_calc_computation_start;
  double vec_calc_computation_end;
  double vec_calc_start = MPI_Wtime();
#endif

  if (rankn == 0) {
    
    int i;
    if (graph_size < comm_size) {

      // Send all nodes if comm size is bigger
      for (i = 0; i < graph_size; i++) {
	MPI_Send(&i, 1, MPI_INT, i+1, tag, MPI_COMM_WORLD);
      }

      // Send termination to finish processes
      int flag;
      for (i = 1; i < comm_size; i++) {
	flag = -1;
        MPI_Send(&flag, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
      }

    } else { // There are more nodes in graph than there are MPI processes

      // First get each process busy
      int send_node; // Corresponds to index of node to send
      int rcv_node;  // Corresponds to name of graph node recieved from worker rank
      vector<int> rcv_gdv;
      for (i = 1; i < comm_size; i++) {

	send_node = (i-1) * CHUNK_SIZE;
        #ifdef DEBUG
          cout << "Sending node " << graph[send_node].Row << " to rank " << i << endl;
        #endif
	MPI_Send(&send_node, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
      }

      // Start probing and recieving results from other processes
      int rec_count = 0;
      send_node += CHUNK_SIZE;
      //int next_job = comm_size-1;
      do {
	
	// First probe for completed work
	int flag;
	MPI_Status master_status;
	MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &master_status);

	if (flag == 1) {
	
	  // Recieve gdv from finished process
	  i = master_status.MPI_SOURCE;
	  int* gdv_array = (int*)calloc(GDV_LENGTH+1, sizeof(int));
	  MPI_Recv(gdv_array, GDV_LENGTH + 1, MPI_INT, i, tag, MPI_COMM_WORLD, &master_status);
	  #ifdef DEBUG
            cout << "Recieved GDV for node " << rcv_node << " from rank " << i << ": " << endl;
	    cout << gdv_array[GDV_LENGTH] << ": ";
            for (int j = 0; j < GDV_LENGTH; j++) {
              cout << gdv_array[j] << ", ";
            }
            cout << endl;
          #endif

	  // Prepare to send next node to finished process.
	  if (send_node < graph_size) { // Jobs still exist.  Send next.
            #ifdef DEBUG
	      cout << "Sending node " << graph[send_node].Row << " to rank " << i << endl;
	    #endif
	    MPI_Send(&send_node, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
	    send_node += 1;
	  } else { // Send termination
	    flag = -1;
	    MPI_Send(&flag, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
	  }

	  // Organize recieved data into returning array
	  rcv_gdv.clear();
          rcv_gdv.resize(GDV_LENGTH);
          for (int j = 0; j < GDV_LENGTH; j++) {
            rcv_gdv[j] = gdv_array[j];
          }
          rcv_node = gdv_array[GDV_LENGTH];

          // We are updating each vertex in the nodes neighborhood.
          // Last GDV update for the node will be node+graph_size 
          // as a signal for when we're done with this node
          if(rcv_node >= graph_size) { // Check if last GDV for the node
            rcv_node -= graph_size;
            rec_count += CHUNK_SIZE;
      	    // Prepare to send next node to finished process.
      	    if (send_node < graph_size) { // Jobs still exist.  Send next.
              #ifdef DEBUG
      	        cout << "Sending node " << graph[send_node].Row << " to rank " << i << endl;
      	      #endif
      	      MPI_Send(&send_node, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
      	      send_node += CHUNK_SIZE;
      	    } else { // Send termination
      	      flag = -1;
      	      MPI_Send(&flag, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
      	    }
          }
//          GDVMetric gdv_extract(rcv_node, rcv_gdv);
//          graph_GDV->push_back(gdv_extract);
//          (*graph_GDV)[rcv_node] = gdv_extract;
          for(int j=0; j<GDV_LENGTH; j++) {
            (*graph_GDV)[rcv_node].GDV[j] += rcv_gdv[j];
          }
      	  free(gdv_array);
      	}
      } while (rec_count < graph_size);

#ifdef RUNTIME
      process_ends_communication = MPI_Wtime();
#endif
      //vec_calc_prior_gather = MPI_Wtime() - vec_calc_start + vec_calc_prior_gather;

      // Sort return vector
      sort(graph_GDV->begin(), graph_GDV->end(), node_id_order());
      #ifdef DEBUG
        cout << "Constructed return GDV array" << endl;
	for (i = 0; i < graph_GDV->size(); i++) {
	  cout << graph_GDV->at(i).node << ": ";
	  for (int j = 0; j < graph_GDV->at(i).GDV.size(); j++) {
	    cout << graph_GDV->at(i).GDV[j] << ", ";
	  }
	  cout << endl;
	}	
      #endif

    }

  }
  else { // Instructions for work processes
    
    int node_name;
    MPI_Status worker_status;

    do {

      MPI_Recv(&node_name, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &worker_status);
      #ifdef DEBUG
	cout << "Recieved node " << graph[node_name].Row << " at rank " << rankn << endl;
      #endif
      if (node_name == -1) {

#ifdef RUNTIME
        process_ends_communication = MPI_Wtime();
#endif
        //vec_calc_prior_gather = MPI_Wtime() - vec_calc_start + vec_calc_prior_gather;
        break;
      }

      int end = node_name+CHUNK_SIZE;
      if(end > graph.size())
        end = graph.size();
      vector<GDVMetric> metrics;
      for(int idx=0; idx<graph.size(); idx++) {
        metrics.push_back(GDVMetric(graph[idx].Row, vector<int>(GDV_LENGTH, 0)));
      }
      for(int node=node_name; node<end; node++) {

#ifdef RUNTIME
        vec_calc_computation_start = MPI_Wtime();
#endif
        Calculate_GDV(graph[node].Row, graph, orbits, metrics);

#ifdef RUTNIME
        if (graph_counter == 1) {
          vec_calc_computation_time_X[node] = MPI_Wtime() - vec_calc_computation_start;
          vec_calc_proc_assign_X[node] = rankn;
        }
        if (graph_counter == 2) {
          vec_calc_computation_time_Y[node] = MPI_Wtime() - vec_calc_computation_start;
          vec_calc_proc_assign_Y[node] = rankn;
        }

//        for(int idx=0; idx<graph.size(); idx++) {
//          int gdv_Length = metrics[idx].GDV.size();
//          int* gdv_array = new int[gdv_Length + 1];
//          for (int j = 0; j < gdv_Length; j++) {
//            gdv_array[j] = metrics[idx].GDV[j];
//          }
//          // Send row+graph_size as a signal that this is the last GDV to update
//          if(idx >= graph.size()-1 && node == end-1) {
//            gdv_array[gdv_Length] = graph[idx].Row + graph.size();
//          } else {
//            gdv_array[gdv_Length] = graph[idx].Row;
//          }
//          MPI_Send(gdv_array, gdv_Length+1, MPI_INT, 0, tag, MPI_COMM_WORLD);
//          delete[] gdv_array;
//        }
#endif

        #ifdef DEBUG
          cout << "Sending GDV for node " << graph[node_name].Row << " from rank " << rankn << endl;
        #endif
      }
      for(int idx=0; idx<graph.size(); idx++) {
        int gdv_Length = metrics[idx].GDV.size();
        int* gdv_array = new int[gdv_Length + 1];
        for (int j = 0; j < gdv_Length; j++) {
          gdv_array[j] = metrics[idx].GDV[j];
        }
        // Send row+graph_size as a signal that this is the last GDV to update
        if(idx >= graph.size()-1) {
          gdv_array[gdv_Length] = graph[idx].Row + graph.size();
        } else {
          gdv_array[gdv_Length] = graph[idx].Row;
        }
        MPI_Send(gdv_array, gdv_Length+1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        delete[] gdv_array;
      }
    } while (node_name != -1); // Exit loop if kill value is sent
/*<<<<<<< HEAD

=======
    graph_GDV->clear();
  }

  //vec_calc_post_gather = MPI_Wtime() - vec_calc_start + vec_calc_post_gather;
#ifdef RUNTIME
  vec_calc_communication_time[graph_counter - 1] = process_ends_communication - vec_calc_start;
#endif
//  cout << "Communication time on process " << rankn << " for graph " << graph_counter << " = " << vec_calc_communication_time[graph_counter - 1] << endl;
  //vec_calc_computation_time = vec_calc_computation_end - vec_calc_computation_start + vec_calc_computation_time;

  #ifdef DEBUG
    cout << "Finished GDV Vector Calc on Rank: " << rankn << endl;
  #endif

}
=======
>>>>>>> 7387cdd... Add framework for passing cmake variables to main

void 
kokkos_GDV_vector_calculation(const matrix_type& graph, 
                              GDVs& graph_GDV, 
                              const Orbits& orbits, 
                              const char* graph_name, 
                              int graph_counter) {
  // Set up parallelization               
  int comm_size, rankn;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankn);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  int tag = 11;
  int graph_size = graph.numRows();
cout << "Starting GDV vector calculation\n";
cout << "# of processes : " << comm_size << endl;
//  int graph_size = graph.size();
//  graph_GDV->clear();
//
//  for(int i=0; i<graph.size(); i++) {
//    vector<int> gdv_vec(GDV_LENGTH, 0);
//    graph_GDV->push_back(GDVMetric(graph[i].Row, gdv_vec));
//  }
  //graph_counter += 1;
  //cout << "Graph counter on rank " << rankn << " = " << graph_counter << endl;

#if RUNTIME
  double process_ends_communication;
  double vec_calc_computation_start;
  double vec_calc_computation_end;
  double vec_calc_start = MPI_Wtime();
#endif

  Kokkos::View<int*> num_combinations("Number of combinations", graph.numRows());
  int k_interval = 1000000;

  Kokkos::View<int*> neighbor_scratch("Neighbors", graph.numRows());
  Kokkos::parallel_for(graph.numRows(), KOKKOS_LAMBDA(const int node) {
    int num_neighbors = 0;
    num_neighbors = EssensKokkos::get_num_neighbors(graph, node, 4, neighbor_scratch);
    for(int i=1; i<5; i++) {
      num_combinations(node) += get_num_combinations(num_neighbors, i);
    }
  });
  printf("Computed # of combinations\n");
  Kokkos::parallel_scan(num_combinations.extent(0), KOKKOS_LAMBDA(const int i, int& update, const bool final) {
    const int val_i = num_combinations(i);
    if(final) {
      num_combinations(i) = update;
    }
    update += val_i;
  });
  Kokkos::View<int*> starts("Start indices", (num_combinations(graph.numRows()-1)/k_interval)+1);
  Kokkos::View<int*> ends("Start indices", (num_combinations(graph.numRows()-1)/k_interval)+1);
  starts(0) = 0;
  int threshold = k_interval;
  int start_index = 1;
  for(int i=0; i<num_combinations.extent(0); i++) {
    if(num_combinations(i) > threshold) {
      starts(start_index++) = i;
      threshold += k_interval;
    }
  }
  for(int i=0; i<starts.extent(0)-1; i++) {
    ends(i) = starts(i+1);
  }
  ends(ends.extent(0)-1) = graph.numRows();
  for(int i=0; i<starts.extent(0); i++) {
    printf("%d ", starts(i));
  }
  printf("\n");
  for(int i=0; i<starts.extent(0); i++) {
    printf("%d ", ends(i));
  }
  printf("\n");
  for(int i=0; i<starts.extent(0); i++) {
    printf("%d ", num_combinations(starts(i)));
  }
  printf("\n");

  if(comm_size == 1)
  {
//    Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>> policy(1, Kokkos::AUTO());
    Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>> policy(1, NUM_THREADS);
    using member_type = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>>::member_type;

cout << "Team size: " << policy.team_size() << endl;
    Kokkos::View<int**> all_neighbors("Neighbor scratch", policy.team_size(), graph.numRows());
    Kokkos::View<int*> combination_counter("Per thread counter", policy.team_size());
    Kokkos::deep_copy(combination_counter, 0);
    Kokkos::Experimental::ScatterView<int**> metrics_sa(graph_GDV);

    Kokkos::View<int**> indices("Index", policy.team_size(), 5);
    Kokkos::View<int**> combination_view("combination", policy.team_size(), 5);
    Kokkos::View<int**> sgraph_distance_signature("dist sig", policy.team_size(), orbits.distance.extent(0));
    Kokkos::View<int**> sgraph_degree_signature("Degree signature", policy.team_size(), 5);
    Kokkos::View<int*[5][5]> induced_subgraph("Subgraph", policy.team_size());
    Kokkos::View<bool** > visited("BFS visited", policy.team_size(), 5);
    Kokkos::View<int**> queue("BFS queue", policy.team_size(), 5);
    Kokkos::View<int**> distance("BFS distance", policy.team_size(), 5);

    int i=0;
    for(i; i<starts.extent(0); i++) {
      Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Dynamic>> bundle_policy(1, NUM_THREADS);
      Kokkos::parallel_for("Calcualte GDV bundle", bundle_policy, KOKKOS_LAMBDA(member_type team_member) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, ends(i)-starts(i)), [=] (int n_offset) {
          int node = starts(i) + n_offset;
          auto neighbor_subview = Kokkos::subview(all_neighbors, team_member.team_rank(), Kokkos::ALL());
          auto indices_subview = Kokkos::subview(indices, team_member.team_rank(), Kokkos::ALL());
          auto combination_subview = Kokkos::subview(combination_view, team_member.team_rank(), Kokkos::ALL());
          auto sgraph_dist_subview = Kokkos::subview(sgraph_distance_signature, team_member.team_rank(), Kokkos::ALL());
          auto sgraph_deg_subview = Kokkos::subview(sgraph_degree_signature, team_member.team_rank(), Kokkos::ALL());
          auto subgraph_subview = Kokkos::subview(induced_subgraph, team_member.team_rank(), Kokkos::ALL(), Kokkos::ALL());
          auto visited_subview = Kokkos::subview(visited, team_member.team_rank(), Kokkos::ALL());
          auto queue_subview = Kokkos::subview(queue, team_member.team_rank(), Kokkos::ALL());
          auto distance_subview = Kokkos::subview(distance, team_member.team_rank(), Kokkos::ALL());
chrono::  steady_clock::time_point t1 = chrono::steady_clock::now();
          kokkos_calculate_GDV(team_member, node, graph, orbits, neighbor_subview, indices_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, combination_counter, 
                              metrics_sa);
          chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
          chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2-t1);
          printf("Done with node: %d, time: %f\n", node, time_span.count());
        });
      });
      Kokkos::Experimental::contribute(graph_GDV, metrics_sa);
      metrics_sa.reset();
    }
//    Kokkos::parallel_for("Calculate GDV", policy, KOKKOS_LAMBDA(member_type team_member) {
//      Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, graph.numRows()), [=] (int node) {
////        int node = team_member.league_rank()*team_member.team_size() + team_member.team_rank();
//        auto neighbor_subview = Kokkos::subview(all_neighbors, team_member.team_rank(), Kokkos::ALL());
//        auto indices_subview = Kokkos::subview(indices, team_member.team_rank(), Kokkos::ALL());
//        auto combination_subview = Kokkos::subview(combination_view, team_member.team_rank(), Kokkos::ALL());
//        auto sgraph_dist_subview = Kokkos::subview(sgraph_distance_signature, team_member.team_rank(), Kokkos::ALL());
//        auto sgraph_deg_subview = Kokkos::subview(sgraph_degree_signature, team_member.team_rank(), Kokkos::ALL());
//        auto subgraph_subview = Kokkos::subview(induced_subgraph, team_member.team_rank(), Kokkos::ALL(), Kokkos::ALL());
//        auto visited_subview = Kokkos::subview(visited, team_member.team_rank(), Kokkos::ALL());
//        auto queue_subview = Kokkos::subview(queue, team_member.team_rank(), Kokkos::ALL());
//        auto distance_subview = Kokkos::subview(distance, team_member.team_rank(), Kokkos::ALL());
//chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
//        kokkos_calculate_GDV(team_member, node, graph, orbits, neighbor_subview, indices_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, combination_counter, 
//                              metrics_sa);
//chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
//chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2-t1);
//        printf("Done with node: %d, time: %f\n", node, time_span.count());
//      });
//    });
//    Kokkos::Experimental::contribute(graph_GDV, metrics_sa);
  }
  else if (rankn == 0) 
  {
    int i;
    if (graph_size < comm_size) {

      // Send all nodes if comm size is bigger
      for (i = 0; i < graph_size; i++) {
        MPI_Send(&i, 1, MPI_INT, i+1, tag, MPI_COMM_WORLD);
      }

      // Send termination to finish processes
      int flag;
      for (i = 1; i < comm_size; i++) {
        flag = -1;
        MPI_Send(&flag, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
      }

    } else { // There are more nodes in graph than there are MPI processes

      // First get each process busy
      int send_node; // Corresponds to index of node to send
      int rcv_node;  // Corresponds to name of graph node recieved from worker rank
      for (i = 1; i < comm_size; i++) {
	      send_node = (i-1) * CHUNK_SIZE;
        #ifdef DEBUG
          cout << "Sending node " << send_node << " to rank " << i << endl;
        #endif
        MPI_Send(&send_node, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
      }
//cout << "Sent initial batch of nodes\n";

      // Start probing and recieving results from other processes
      int rec_count = 0;
      send_node += CHUNK_SIZE;
      //int next_job = comm_size-1;
      do {
	
        // First probe for completed work
      	int flag;
      	MPI_Status master_status;
      	MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &master_status);
      
      	if (flag == 1) {
      	
      	  // Recieve gdv from finished process
      	  i = master_status.MPI_SOURCE;
          Kokkos::View<int[GDV_LENGTH+1]> gdv_array_d("Recv buffer");
          auto gdv_array_h = Kokkos::create_mirror_view(gdv_array_d);
      	  MPI_Recv(gdv_array_h.data(), GDV_LENGTH + 1, MPI_INT, i, tag, MPI_COMM_WORLD, &master_status);
          Kokkos::deep_copy(gdv_array_d, gdv_array_h);
//cout << "Received GDV from " << gdv_array_h(GDV_LENGTH) << endl;
      	  #ifdef DEBUG
            cout << "Recieved GDV for node " << rcv_node << " from rank " << i << ": " << endl;
            cout << gdv_array_h(GDV_LENGTH) << ": ";
            for (int j = 0; j < GDV_LENGTH; j++) {
              cout << gdv_array_h(j) << ", ";
            }
            cout << endl;
          #endif
//cout << "Receive buffer for node " << gdv_array_h(GDV_LENGTH) << ": ";
//for(int k=0; k<gdv_array_h.extent(0); k++) {
//  cout << gdv_array_h(k) << " ";
//}
//cout << endl;
      
      	  // Organize recieved data into returning array
//      	  rcv_gdv.clear();
//          rcv_gdv.resize(GDV_LENGTH);
//          for (int j = 0; j < GDV_LENGTH; j++) {
//            rcv_gdv[j] = gdv_array[j];
//          }
//          Kokkos::parallel_for("Copy GDV", Kokkos::RangePolicy<>(0, GDV_LENGTH), KOKKOS_LAMBDA(const int j) {
//            rcv_gdv(j) = gdv_array_d[j];
//          });
          rcv_node = gdv_array_h(GDV_LENGTH);
          // We are updating each vertex in the nodes neighborhood.
          // Last GDV update for the node will be node+graph_size 
          // as a signal for when we're done with this node
          if(rcv_node >= graph_size) { // Check if last GDV for the node
            rcv_node -= graph_size;
            rec_count += CHUNK_SIZE;
//cout << "Received GDVs for node " << rcv_node << endl;
      	    // Prepare to send next node to finished process.
      	    if (send_node < graph_size) { // Jobs still exist.  Send next.
              #ifdef DEBUG
//      	        cout << "Sending node " << graph[send_node].Row << " to rank " << i << endl;
      	        cout << "Sending node " << send_node << " to rank " << i << endl;
      	      #endif
      	      MPI_Send(&send_node, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
      	      send_node += CHUNK_SIZE;
      	    } else { // Send termination
      	      flag = -1;
      	      MPI_Send(&flag, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
      	    }
//cout << "Send new batch of nodes to " << i << endl;
          }
//          for(int j=0; j<GDV_LENGTH; j++) {
//            (*graph_GDV)[rcv_node].GDV[j] += rcv_gdv[j];
//          }
          for(int j=0; j<GDV_LENGTH; j++) {
            graph_GDV(rcv_node, j) += gdv_array_d(j);
          }
//cout << "Updated GDV for " << rcv_node << endl;
//
//cout << "GDV: Root node " << rcv_node << endl;
//for(int i=0; i<graph_GDV.extent(0); i++) {
//  cout << "Node " << i << ": ";
//  for(int j=0; j<graph_GDV.extent(1); j++) {
//    cout << graph_GDV(i,j) << " ";
//  }
//  cout << endl;
//}
//cout << endl;
//      	  free(gdv_array);
      	}
      } while (rec_count < graph_size);

#if RUNTIME
      process_ends_communication = MPI_Wtime();
      //vec_calc_prior_gather = MPI_Wtime() - vec_calc_start + vec_calc_prior_gather;
#endif

      // Sort return vector
//      sort(graph_GDV->begin(), graph_GDV->end(), node_id_order());
      #ifdef DEBUG
        cout << "Constructed return GDV array" << endl;
        for (i = 0; i < graph_GDV->size(); i++) {
      	  cout << graph_GDV->at(i).node << ": ";
      	  for (int j = 0; j < graph_GDV->at(i).GDV.size(); j++) {
      	    cout << graph_GDV->at(i).GDV[j] << ", ";
      	  }
      	  cout << endl;
      	}	
      #endif
    }
  }
  else // Instructions for work processes
  { 
    int node_name;
    MPI_Status worker_status;

    do {

//cout << "Starting worker thread\n";
      MPI_Recv(&node_name, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &worker_status);
//cout << "Got new batch of nodes starting at " << node_name << "\n";
      #ifdef DEBUG
        cout << "Recieved node " << node_name << " at rank " << rankn << endl;
      #endif
      if (node_name == -1) {
#if RUNTIME
        process_ends_communication = MPI_Wtime();
#endif
        //vec_calc_prior_gather = MPI_Wtime() - vec_calc_start + vec_calc_prior_gather;
//cout << "Thread " << rankn << " done\n";
        break;
      }
      int end = node_name+CHUNK_SIZE;
//      if(end > graph.size())
//        end = graph.size();
      int nrows = graph.numRows();
      if(end > nrows)
        end = nrows;
//      vector<GDVMetric> metrics;
//      for(int idx=0; idx<graph.numRows(); idx++) {
//        metrics.push_back(GDVMetric(graph[idx].Row, vector<int>(GDV_LENGTH, 0)));
//      }

//      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(1, Kokkos::AUTO());
//      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(end-node_name, Kokkos::AUTO());
      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(1, NUM_THREADS);
//      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy((end-node_name)/4, 4);
//      policy.set_scratch_size(0, Kokkos::PerThread(128));
      using member_type = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;

cout << "Team size: " << policy.team_size() << endl;
      Kokkos::View<int**> all_neighbors("Neighbor scratch", policy.team_size(), graph.numRows());
      GDVs metrics("GDVs", graph.numRows(), GDV_LENGTH);
      Kokkos::Experimental::ScatterView<int**> metrics_sa(metrics);

  Kokkos::View<int*> combination_counter("Per thread counter", policy.team_size());
  Kokkos::deep_copy(combination_counter, 0);
  Kokkos::View<int** > indices("Index", policy.team_size(), 5);
  Kokkos::View<int** > combination_view("combination", policy.team_size(), 5);
  Kokkos::View<int** > sgraph_distance_signature("dist sig", policy.team_size(), orbits.distance.extent(0));
  Kokkos::View<int**> sgraph_degree_signature("Degree signature", policy.team_size(), 5);
  Kokkos::View<int*[5][5]> induced_subgraph("Subgraph", policy.team_size());
  Kokkos::View<bool** > visited("BFS visited", policy.team_size(), 5);
  Kokkos::View<int**> queue("BFS queue", policy.team_size(), 5);
  Kokkos::View<int**> distance("BFS distance", policy.team_size(), 5);
//      Kokkos::View<int[1]> node_counter("Node counter");
//      node_counter(0) = node_name;
      Kokkos::parallel_for("Calculate GDV", policy, KOKKOS_LAMBDA(member_type team_member) {
//        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, end-node_name), [=] (int n) {
          int node = node_name + team_member.league_rank()*team_member.team_size() + team_member.team_rank();
//          int node = Kokkos::atomic_fetch_add(&node_counter(0), 1);
//          int node = node_name + n;
//printf("Node: %d\tThreadid: %d\tTeam size: %d\tLeague size: %d\n", node, team_member.team_rank(), team_member.team_size(), team_member.league_size());
          auto neighbor_subview = Kokkos::subview(all_neighbors, team_member.team_rank(), Kokkos::ALL());
          auto indices_subview = Kokkos::subview(indices, team_member.team_rank(), Kokkos::ALL());
          auto combination_subview = Kokkos::subview(combination_view, team_member.team_rank(), Kokkos::ALL());
          auto sgraph_dist_subview = Kokkos::subview(sgraph_distance_signature, team_member.team_rank(), Kokkos::ALL());
          auto sgraph_deg_subview = Kokkos::subview(sgraph_degree_signature, team_member.team_rank(), Kokkos::ALL());
          auto subgraph_subview = Kokkos::subview(induced_subgraph, team_member.team_rank(), Kokkos::ALL(), Kokkos::ALL());
          auto visited_subview = Kokkos::subview(visited, team_member.team_rank(), Kokkos::ALL());
          auto queue_subview = Kokkos::subview(queue, team_member.team_rank(), Kokkos::ALL());
          auto distance_subview = Kokkos::subview(distance, team_member.team_rank(), Kokkos::ALL());
          kokkos_calculate_GDV(team_member, node, graph, orbits, neighbor_subview, indices_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, combination_counter, 
                                metrics_sa);
//        });
      });
      Kokkos::Experimental::contribute(metrics, metrics_sa);
      Kokkos::View<int*> gdv_array("Send buffer for GDV", GDV_LENGTH+1);
      for(int idx=0; idx<nrows; idx++) {
        for(int j=0; j<GDV_LENGTH; j++) {
          gdv_array(j) = metrics(idx, j);
        }
        if(idx >= nrows-1) {
          gdv_array(GDV_LENGTH) = idx + nrows;
        } else {
          gdv_array(GDV_LENGTH) = idx;
        }
//cout << "Send buffer for node " << idx << ": ";
//for(int k=0; k<gdv_array.extent(0); k++) {
//  cout << gdv_array(k) << " ";
//}
//cout << endl;
        MPI_Send(gdv_array.data(), GDV_LENGTH+1, MPI_INT, 0, tag, MPI_COMM_WORLD);
      }
    } while (node_name != -1); // Exit loop if kill value is sent
//    graph_GDV->clear();
>>>>>>> f86ac41... Updating runtime tracking to only happen on ifdef*/
  }

  //vec_calc_post_gather = MPI_Wtime() - vec_calc_start + vec_calc_post_gather;
#if RUNTIME
  vec_calc_communication_time[graph_counter - 1] = process_ends_communication - vec_calc_start;
#endif
//  cout << "Communication time on process " << rankn << " for graph " << graph_counter << " = " << vec_calc_communication_time[graph_counter - 1] << endl;
  //vec_calc_computation_time = vec_calc_computation_end - vec_calc_computation_start + vec_calc_computation_time;

  #ifdef DEBUG
    cout << "Finished GDV Vector Calc on Rank: " << rankn << endl;
  #endif

}


void Calculate_GDV(int node,A_Network Graph,vector<OrbitMetric> &orbits, vector<GDVMetric> &gdvMetrics)
{
  GDV_functions gdvf;
//  vector<int> gdv(orbits.size(),0);
  // printf("calculating GDV for node %d\n",node);
  vector<int> neighbours;
  gdvf.find_neighbours(node,Graph,4,&neighbours);
//if(node == 100) {
//  printf("Neighbours of node 100: ");
//  for(int i=0; i<neighbours.size(); i++) {
//    printf("%d ", neighbours[i]);
//  }
//  printf("\n");
//}
  //print_vector(neighbours);
  int set[neighbours.size()]; 
  std::copy( neighbours.begin(), neighbours.end(), set );
  //int numElements = *(&set + 1) - set;
  int numElements = sizeof(set)/sizeof(set[0]);
  for (int node_count = 1; node_count < 5; node_count++)
  {
//    vector<vector<int>> combinationsList;
//    gdvf.find_combinations(set, numElements,node_count,&combinationsList);
    CombinationGenerator generator(neighbours.size(), node_count);
    vector<int> combination(node_count, 0);
//printf("Size of combination: %d\n", combination.size());
//printf("Number of neighbours: %d\n", neighbours.size());
//printf("Neighbor 0: %d\n", set[0]);
//    for (vector<int> &combination : combinationsList)
//    {
    while(!generator.done) {
      generator.get_combination(neighbours, combination);
      A_Network induced_sgraph;
      vector<int> subgraph_degree_signature;
      vector<int> subgraph_distance_signature;
      bool is_connected = false;
      combination.push_back(node);
//printf("Combination: ");
//for(int i=0; i<combination.size(); i++) {
//  printf("%d ", combination[i]);
//}
//printf("\n");
      gdvf.inducedSubgraph(Graph, combination, induced_sgraph);
      gdvf.isConnected(induced_sgraph, is_connected);
      if(is_connected)
      {
        gdvf.degree_signature(induced_sgraph,subgraph_degree_signature);
        sort(subgraph_degree_signature.begin(), subgraph_degree_signature.end());
        vector<OrbitMetric> filter_orbits;
        gdvf.orbit_filter(orbits, node_count+1, filter_orbits);
        for(int idx=0; idx<induced_sgraph.size(); idx++)
        {
          int v = induced_sgraph[idx].Row;
          gdvf.distance_signature(v,induced_sgraph,subgraph_distance_signature);
          for(OrbitMetric orbit: filter_orbits)
          {
            if( orbit.orbitDistance == subgraph_distance_signature && 
                orbit.orbitDegree == subgraph_degree_signature)
            {
              gdvMetrics[v].GDV[orbit.orbitNumber] += 1;
              break;
            }
          }
        }
      }
      generator.next();
    }
  }
}

/*void Calculate_GDV(int node,A_Network Graph,vector<OrbitMetric> &orbits, GDVMetric &gdvMetric)
{
    GDV_functions gdvf;
    vector<int> gdv(orbits.size(),0);
    // printf("calculating GDV for node %d\n",node);
    vector<int> neighbours;
    gdvf.find_neighbours(node,Graph,4,&neighbours);
    //print_vector(neighbours);
    int set[neighbours.size()]; 
    std::copy( neighbours.begin(), neighbours.end(), set );
    //int numElements = *(&set + 1) - set;
    int numElements = sizeof(set)/sizeof(set[0]);
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
}*/

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


/*
bool sort_by_dest(vector<int>& a, vector<int>& b) {
  return a[1] < b[1];
}

KOKKOS_INLINE_FUNCTION void readin_graph(ifstream* file, matrix_type& graph) 
{
  string line;;
  int nnz = 0;

  int lastrow = 0;
  vector<int> rows, cols;
  vector<float> vals;
  int k = 0;
  vector<vector<int>> edge_list;
  while(std::getline(*file,line))
  {
    int u, v;
    float w;
    sscanf(line.c_str(), "%d %d %f", &u, &v, &w);
    //printf("Found nodes u and v: %d, %d\n", u, v);
    if(u > lastrow) 
      lastrow = u;
    if(v > lastrow) 
      lastrow = v;
    vector<int> edge1;
    edge1.push_back(u);
    edge1.push_back(v);
    if(w != 0.0) {
        edge1.push_back(1);
    } else {
        edge1.push_back(0);
    }
    edge_list.push_back(edge1);
    if(u != v) {
      vector<int> edge2;
      edge2.push_back(v);
      edge2.push_back(u);
      if(w != 0.0) {
            edge2.push_back(1);
      } else {
            edge2.push_back(0);
      }
      edge_list.push_back(edge2);
    }
    //printf("Edge list size is: %d\n", edge_list.size());
    //printf("Last row is: %d\n", lastrow); 
  }
  sort(edge_list.begin(), edge_list.end(), sort_by_leading_edge);

  vector<int> rowmap(lastrow+2, 0);
  for(int i=0; i<edge_list.size(); i++) {
    rows.push_back(edge_list[i][0]);
    cols.push_back(edge_list[i][1]);
    vals.push_back(edge_list[i][2]);
  }
  
  for(int i=0; i<rows.size(); i++) {
    if(rows[i] != cols[i])
    {
      rowmap[rows[i]+1]++;
    }
    else
    {
      rowmap[rows[i]+1]++;
    }
  }
  for(int i=1; i<rowmap.size(); i++) {
    rowmap[i] += rowmap[i-1];
  }
  for(int i=0; i<rowmap.size()-1; i++) {
    sort(cols.begin()+rowmap[i], cols.begin()+rowmap[i+1]);
  }
  
  graph = matrix_type("Graph", lastrow+1, lastrow+1, vals.size(), 
                      vals.data(), rowmap.data(), cols.data());
  return;
}
*/




