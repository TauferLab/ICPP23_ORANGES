// Required directives are to be added from ESSENS to use the functionality 

#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
#include "../headers/GDV_functions.hpp"
#include "../headers/class_definitions.hpp"
#include "../headers/combinations.hpp"
#include "../headers/kokkos_functions.hpp"
#include "../headers/dirty_page_tracking.hpp"
#include "../headers/hash_change_detection.hpp"
#include <time.h>
#include <stdlib.h>
#include <ctime>
#include <mpi.h>
#include <fstream>
#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include <Kokkos_ScatterView.hpp>
#include <Kokkos_UnorderedMap.hpp>

#define GDV_LENGTH 73
#define CHUNK_SIZE 1
//#define NUM_THREADS Kokkos::AUTO
#define NUM_THREADS 4
#define MAX_SUBGRAPH 5
//#define ANALYSIS
//#define HASH_ANALYSIS
//#define SCAN_ANALYSIS
//#define DEBUG
//#define DIRTY_PAGE_TRACKING
//#define RESILIENCE
//#define AUTO_CHECKPOINT
//#define OUTPUT_MATRIX
//#define OUTPUT_GDV
//#define OUTPUT_MATRIX
#define LAYOUT Kokkos::LayoutRight

#ifdef AUTO_CHECKPOINT
#include <resilience/Resilience.hpp>
#include <resilience/CheckpointFilter.hpp>
#include <veloc.h>
#endif

void Calculate_GDV(int ,A_Network ,vector<OrbitMetric>&, GDVMetric&);
void Calculate_GDV(int node,A_Network Graph,vector<OrbitMetric> &orbits, vector<GDVMetric> &gdvMetrics);
void readin_orbits(  ifstream* ,vector<OrbitMetric>* );
void convert_string_vector_int(string* , vector<int>* ,string );
void Similarity_Metric_calculation_for_two_graphs(A_Network, A_Network,vector<OrbitMetric>, string, string);
double GDV_distance_calculation(GDVMetric&, GDVMetric&);
void metric_formula(GDVMetric&, double*);
void GDV_vector_calculation(A_Network,vector<GDVMetric>*,  vector<OrbitMetric>, const char*, int);

void kokkos_readin_orbits(ifstream *file, Orbits& orbits );
void readin_graph(ifstream* file, matrix_type& graph);
void kokkos_Similarity_Metric_calculation_for_two_graphs(const matrix_type&, const matrix_type&, const Orbits&, string, string);
template <class SubviewType>
KOKKOS_INLINE_FUNCTION double kokkos_GDV_distance_calculation(SubviewType&, SubviewType&);
template <class SubviewType>
KOKKOS_INLINE_FUNCTION double kokkos_metric_formula(SubviewType &gdvm);
void kokkos_GDV_vector_calculation(const matrix_type&, GDVs&, const Orbits&, const char*, int);
template<class NeighborView, class IntView, class GraphView, class BoolView, class Scatter>
KOKKOS_INLINE_FUNCTION void 
kokkos_calculate_GDV(Kokkos::TeamPolicy<>::member_type team_member,
                      int node, 
                      const matrix_type& graph, 
                      const Orbits& orbits, 
                      NeighborView& neighbor_buff,
                      IntView& indices,
                      IntView& combination_view,
                      IntView& sgraph_distance_signature,
                      IntView& sgraph_degree_signature,
                      GraphView& induced_subgraph,
                      BoolView& visited,
                      IntView& queue,
                      IntView& distance,
                      Scatter gdvMetrics_sa
                    );
template<class NeighborView, class IntView, class GraphView, class BoolView, class Scatter>
KOKKOS_INLINE_FUNCTION void 
kokkos_calculate_GDV(Kokkos::TeamPolicy<>::member_type team_member,
                      int node, 
                      int node_count,
                      const matrix_type& graph, 
                      const Orbits& orbits, 
                      const NeighborView& neighbors,
                      const int n_neighbors,
                      const IntView& combination_view,
                      IntView& sgraph_distance_signature,
                      IntView& sgraph_degree_signature,
                      GraphView& induced_subgraph,
                      BoolView& visited,
                      IntView& queue,
                      IntView& distance,
                      Scatter gdvMetrics_sa
                    );

struct node_id_order {
  inline bool operator() (const GDVMetric& gdv_met_1, const GDVMetric& gdv_met_2) {
    return (gdv_met_1.node < gdv_met_2.node);
  }
};

// Define variables for keeping track of time for load imbalancing tests.
double total_time_taken;
double vec_calc_communication_time[2] = {};
double pre_process_time;
double report_output_time;
double* vec_calc_computation_time_X;
double* vec_calc_computation_time_Y;
int* vec_calc_proc_assign_X;
int* vec_calc_proc_assign_Y;
//double vec_calc_prior_gather;
//double vec_calc_post_gather;
int vec_calc_avg_node_deg;

using namespace std;

int main(int argc, char *argv[]) {

//  MPI_Init(&argc,&argv);
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided); 
  Kokkos::initialize(argc, argv);
  {
  //  if (rank == 0) {
  //clock_t out_tStart = clock();
  time_t now = time(0);
    //}
  clock_t tStart = clock();
  clock_t q, q1, q2,t;
  double pre_tStart = MPI_Wtime();
//  GDV_functions gdvf;
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

  Orbits k_orbits;
  kokkos_readin_orbits(&the_file2, k_orbits);

  matrix_type graphx, graphy;
  readin_graph(&the_file0, graphx);
  readin_graph(&the_file1, graphy);

  //clock_t out_tStart = clock();

  int numtasks, rank, dest, source, rc, count, tag=0;
  MPI_Status Stat;   // required variable for receive routines                                                                                                                                                                          

  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) 
    printf("Comm size: %d\n", numtasks);

  // Allocate space for time recording by graph node
  vec_calc_computation_time_X = new double[graphx.numRows()]();
  vec_calc_computation_time_Y = new double[graphy.numRows()]();
  vec_calc_proc_assign_X = new int[graphx.numRows()]();
  vec_calc_proc_assign_Y = new int[graphy.numRows()]();

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
  if(rank == 0) {
    cout << "# of vertices: " << graphx.numRows() << endl;
    cout << "# of vertices: " << graphy.numRows() << endl;
  }

  pre_process_time = (double)(MPI_Wtime() - pre_tStart);

  // Perform Similarity Calculations
//  Similarity_Metric_calculation_for_two_graphs(X,Y,orbits, graph_name1, graph_name2);
  kokkos_Similarity_Metric_calculation_for_two_graphs(graphx, graphy, k_orbits, graph_name1, graph_name2);

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

Kokkos::TeamPolicy<> team_policy(1, NUM_THREADS);
  
    if (myfile.is_open()) {
      myfile << month << "/" << day << "/" << year << "\t" << numtasks << "\t" << team_policy.team_size() << "\t" << graph_name1 << "\t" << graph_name2 << "\t" << graphx.numRows() << "\t\t" << total_time_taken << " \n";
      myfile.close();
    }
    else { 
      cout << "Out File Did Not Open";
    }

    // Print out rank specific runtime data
    string time_file = "runtime_data/runtimes_rec_over.txt";
    myfile.open(time_file, ofstream::trunc);
    /*if (!myfile.is_open()) {
      cout << "INPUT ERROR:: Could not open the local time recording file\n";
    }
    myfile << "Time Taken in Similarity Metric Calculation = " << " \n";
    myfile << "Rank\tGraph 1\tGraph 2\tTotal\n";
    for (int i = 0; i < numtasks; i++) {
      myfile << i << " " << time_buff[num_times*i+1] << " " << time_buff[num_times*i+2] << " " << time_buff[num_times*i] << " \n";
    }
    myfile.close();*/
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
    for (int i = 0; i < graphx.numRows(); i++) {
      myfile << i << " " << vec_calc_proc_assign_X[i] << " " << vec_calc_computation_time_X[i] << endl;
    }
    myfile.close();
  }
  myfile.open(computation_time_file_y, ofstream::trunc);
  if (!myfile.is_open()) {
    cout << "INPUT ERROR:: Could not open the local time recording file\n";
  }
  if (myfile.is_open()) {
    for (int i = 0; i < graphy.numRows(); i++) {
      myfile << i << " " << vec_calc_proc_assign_Y[i] << " " << vec_calc_computation_time_Y[i] << endl;
    }
    myfile.close();
  }


  //printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  //cout << "Time taken on rank " << rank << " = " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
  //MPI_Barrier( MPI_COMM_WORLD );
  delete[] vec_calc_computation_time_X;
  delete[] vec_calc_computation_time_Y;
  delete[] vec_calc_proc_assign_X;
  delete[] vec_calc_proc_assign_Y;
  }
  Kokkos::finalize();
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

  int graph_counter = 1;
  GDV_vector_calculation(graph1, &graph1_GDV, orbits, "graph1", graph_counter); 

  graph_counter = 2;
  GDV_vector_calculation(graph2, &graph2_GDV, orbits, "graph2", graph_counter); 
  MPI_Barrier(MPI_COMM_WORLD);

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

KOKKOS_INLINE_FUNCTION void 
kokkos_Similarity_Metric_calculation_for_two_graphs(const matrix_type& graph1, 
                                                    const matrix_type& graph2, 
                                                    const Orbits& orbits, 
                                                    string graph_tag1, 
                                                    string graph_tag2) {
  //clock_t out_tStart = clock();
  MPI_Barrier( MPI_COMM_WORLD );
  double out_tStart = MPI_Wtime();

  int num_orbits = orbits.degree.extent(0);
  GDVs graph1_GDV("GDV for graph1", graph1.numRows(), num_orbits);
  GDVs graph2_GDV("GDV for graph2", graph2.numRows(), num_orbits);

  int rankm, numtasksm;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankm);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasksm);

  int graph_counter = 1;
//  kokkos_GDV_vector_calculation(graph1, graph1_GDV, orbits, "graph1", graph_counter); 
  kokkos_GDV_vector_calculation(graph1, graph1_GDV, orbits, graph_tag1.c_str(), graph_counter); 
  MPI_Barrier(MPI_COMM_WORLD);

  if(rankm == 0)
    cout << "Done with GDV calculation for graph 1\n";

#ifdef DEBUG
  if(rankm == 0) {
    cout << "Post kokkos_GDV_vector_calculation: graph 1" << endl;
    for(int i=0; i<graph1_GDV.extent(0); i++) {
      cout << "Node " << i << ": ";
      for(int j=0; j<graph1_GDV.extent(1); j++) {
        cout << graph1_GDV(i,j) << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
#endif

#ifdef AUTO_CHECKPOINT
VELOC_Checkpoint_wait();
#endif

  graph_counter = 2;
//  kokkos_GDV_vector_calculation(graph2, graph2_GDV, orbits, "graph2", graph_counter); 
  kokkos_GDV_vector_calculation(graph2, graph2_GDV, orbits, graph_tag2.c_str(), graph_counter); 
  MPI_Barrier(MPI_COMM_WORLD);

  if(rankm == 0)
    cout << "Done with GDV calculation for graph 2\n";

#ifdef DEBUG
  if(rankm == 0) {
    cout << "Post kokkos_GDV_vector_calculation: graph 2" << endl;
    for(int i=0; i<graph2_GDV.extent(0); i++) {
      cout << "Node " << i << ": ";
      for(int j=0; j<graph2_GDV.extent(1); j++) {
        cout << graph2_GDV(i,j) << " ";
      }
      cout << endl;
    }
    cout << endl;
  }
#endif

#ifdef OUTPUT_MATRIX
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
#endif

  // Measure final total time 
  total_time_taken = (double)(MPI_Wtime() - out_tStart);

  #ifdef DEBUG
    cout << "Finished similarity matrix calculation and prepped for fileIO on Rank: " << rankm << endl;
  #endif

  double report_tStart = MPI_Wtime();

  // Perform File I/O for Output Data
  if (rankm == 0) {

    // Start by Printing out Similarity Matrix into a file
    ofstream myfile; 
#ifdef OUTPUT_MATRIX
    string filename = "out_similarity_matrix.txt";
    myfile.open(filename);
    for(int i=1; i<m;i++) {
      for(int j=1;j<n;j++) {
        myfile<<" { "<<sim_mat(i,j)<<" } ";
      }
      myfile<<"||"<<endl;
    }
    myfile.close();
#endif
    
#ifdef OUTPUT_GDV    
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
#endif
  }
  report_output_time = (double)(MPI_Wtime() - report_tStart);
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

template <class SubviewType>
KOKKOS_INLINE_FUNCTION
double kokkos_GDV_distance_calculation(SubviewType& gdv1, SubviewType& gdv2)
{
  double gdv1_score;
  double gdv2_score;
  gdv1_score = kokkos_metric_formula(gdv1);
  gdv2_score = kokkos_metric_formula(gdv2);
  
  double similarity_score = abs(gdv1_score - gdv2_score);
  return similarity_score;
}

template <class SubviewType>
KOKKOS_INLINE_FUNCTION
double kokkos_metric_formula(SubviewType &gdvm)
{
  double sum = 0;
  for(size_t i=0; i<gdvm.extent(0); i++) {
    sum += static_cast<double>(static_cast<int64_t>(gdvm[i])*static_cast<int64_t>(gdvm[i]));
  }
  return sqrt(sum);
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
  for(int i=0; i<graph.size(); i++) {
    vector<int> gdv_vec(GDV_LENGTH, 0);
    graph_GDV->push_back(GDVMetric(graph[i].Row, gdv_vec));
  }

  double process_ends_communication;
  double vec_calc_computation_start;
  double vec_calc_computation_end;
  double vec_calc_start = MPI_Wtime();

  if (rankn == 0) 
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
      
//      	  // Prepare to send next node to finished process.
//      	  if (send_node < graph_size) { // Jobs still exist.  Send next.
//            #ifdef DEBUG
//      	      cout << "Sending node " << graph[send_node].Row << " to rank " << i << endl;
//      	    #endif
//      	    MPI_Send(&send_node, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
//      	    send_node += 1;
//      	  } else { // Send termination
//      	    flag = -1;
//      	    MPI_Send(&flag, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
//      	  }
      
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
          for(int j=0; j<GDV_LENGTH; j++) {
            (*graph_GDV)[rcv_node].GDV[j] += rcv_gdv[j];
          }
      	  free(gdv_array);
      	}
      } while (rec_count < graph_size);

      process_ends_communication = MPI_Wtime();
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
  else // Instructions for work processes
  { 
    int node_name;
    MPI_Status worker_status;

    do {

      MPI_Recv(&node_name, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &worker_status);
      #ifdef DEBUG
        cout << "Recieved node " << graph[node_name].Row << " at rank " << rankn << endl;
      #endif
      if (node_name == -1) {
        process_ends_communication = MPI_Wtime();
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
        vec_calc_computation_start = MPI_Wtime();
        Calculate_GDV(graph[node].Row, graph, orbits, metrics);

        if (graph_counter == 1) {
          vec_calc_computation_time_X[node] = MPI_Wtime() - vec_calc_computation_start;
          vec_calc_proc_assign_X[node] = rankn;
        }
        if (graph_counter == 2) {
          vec_calc_computation_time_Y[node] = MPI_Wtime() - vec_calc_computation_start;
          vec_calc_proc_assign_Y[node] = rankn;
        }
        #ifdef DEBUG
          cout << "Sending GDV for node " << graph[node_name].Row << " from rank " << rankn << endl;
        #endif
      }
      vector<int> gdv_array(metrics[0].GDV.size()+1, 0);
      for(int idx=0; idx<graph.size(); idx++) {
        int gdv_Length = metrics[idx].GDV.size();
        gdv_array.resize(gdv_Length+1);
        for (int j = 0; j < gdv_Length; j++) {
          gdv_array[j] = metrics[idx].GDV[j];
        }
        // Send row+graph_size as a signal that this is the last GDV to update
        if(idx >= graph.size()-1) {
          gdv_array[gdv_Length] = graph[idx].Row + graph.size();
        } else {
          gdv_array[gdv_Length] = graph[idx].Row;
        }
        MPI_Send(gdv_array.data(), gdv_Length+1, MPI_INT, 0, tag, MPI_COMM_WORLD);
      }
    } while (node_name != -1); // Exit loop if kill value is sent
    graph_GDV->clear();
  }

  //vec_calc_post_gather = MPI_Wtime() - vec_calc_start + vec_calc_post_gather;
  vec_calc_communication_time[graph_counter - 1] = process_ends_communication - vec_calc_start;
//  cout << "Communication time on process " << rankn << " for graph " << graph_counter << " = " << vec_calc_communication_time[graph_counter - 1] << endl;
  //vec_calc_computation_time = vec_calc_computation_end - vec_calc_computation_start + vec_calc_computation_time;

  #ifdef DEBUG
    cout << "Finished GDV Vector Calc on Rank: " << rankn << endl;
  #endif

}

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
  uint32_t graph_size = graph.numRows();
  std::string label(graph_name);
  label += "_Rank";
  label += std::to_string(rankn);

#ifdef DIRTY_PAGE_TRACKING
  pid_t pid = getpid();
  uint64_t page_size = sysconf(_SC_PAGE_SIZE);
#endif

#ifdef DEBUG
  cout << "Starting GDV vector calculation\n";
  cout << "# of processes : " << comm_size << endl;
#endif

  double process_ends_communication;
  double vec_calc_computation_start;
  double vec_calc_computation_end;
  double vec_calc_start = MPI_Wtime();

  Kokkos::View<uint64_t*> num_combinations("Number of combinations", graph.numRows());
  Kokkos::View<uint32_t*> num_neighbors("Number of neighbors", graph.numRows());
//  uint64_t k_interval = 2400000000;
  uint64_t k_interval = 1000000000;
//  uint64_t k_interval = 200000000;
//  uint64_t k_interval = 1000000;

  Kokkos::TeamPolicy<> team_policy(1, NUM_THREADS);
#ifdef DEBUG
  printf("Number of threads: %d\n", team_policy.team_size());
#endif
  using team_member_type = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Dynamic>>::member_type;
  Kokkos::View<int**> neighbor_scratch("Neighbors", team_policy.team_size(), graph.numRows());
#ifdef DEBUG
  printf("Starting combination counter kernel\n");
#endif

  uint32_t vertices_per_proc = graph.numRows()/comm_size;
  size_t start_offset = rankn*vertices_per_proc;
  size_t end_offset = start_offset+vertices_per_proc;
  if(end_offset > graph.numRows())
    end_offset = graph.numRows();
  int* recv_count = new int[comm_size];
  int* displs = new int[comm_size];
  int running_count = 0;
  for(size_t i=0; i<comm_size-1; i++) {
    recv_count[i] = vertices_per_proc;
    displs[i] = running_count;
    running_count += vertices_per_proc;
  }
  recv_count[size_t (comm_size-1)] = graph.numRows()-(comm_size-1)*vertices_per_proc;
  displs[static_cast<size_t>(comm_size-1)] = running_count;
#ifdef DEBUG
  for(int i=0; i<comm_size; i++) {
    printf("(%d,%d) ", recv_count[i], displs[i]);
  }
  printf("\n");
  printf("(%u,%u)\n", start_offset, end_offset);
#endif
  Kokkos::parallel_for("Get # of combinations", team_policy, KOKKOS_LAMBDA(team_member_type team_member) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, end_offset-start_offset), [=] (int n) {
      size_t node = start_offset+n;
      auto neighbor_subview = Kokkos::subview(neighbor_scratch, team_member.team_rank(), Kokkos::ALL());
      uint32_t n_neighbors = 0;
      n_neighbors = EssensKokkos::get_num_neighbors(graph, node, 4, neighbor_subview);
//      n_neighbors = get_approx_num_neighbors(graph, node, 4);
      num_neighbors(node) = n_neighbors;
      for(int i=1; i<5; i++) {
        num_combinations(node) += get_num_combinations(n_neighbors, i);
      }
    });
  });
#ifdef DEBUG
printf("Calculated # of combinations\n");
#endif
  uint64_t total_combinations = 0;
  uint64_t local_combinations = 0;
  Kokkos::parallel_reduce("Number of combinations", end_offset-start_offset, 
  KOKKOS_LAMBDA(const int i, uint64_t& update) {
    update += num_combinations(i+start_offset);
  }, local_combinations);

  Kokkos::View<uint64_t*> send_comb("Send buffer", end_offset-start_offset);
  auto s_view = Kokkos::subview(num_combinations, Kokkos::make_pair(start_offset, end_offset));
  Kokkos::deep_copy(send_comb, s_view);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allgatherv(send_comb.data(), 
                  end_offset-start_offset, 
                  MPI_UNSIGNED_LONG, 
                  num_combinations.data(), 
                  recv_count, 
                  displs, 
                  MPI_UNSIGNED_LONG, 
                  MPI_COMM_WORLD);
  delete[] recv_count, displs;
  MPI_Allreduce(&local_combinations, &total_combinations, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

//  Kokkos::parallel_for("Get # combinations", team_policy, KOKKOS_LAMBDA(team_member_type team_member) {
//    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, graph.numRows()), [=] (int node) {
//      auto neighbor_subview = Kokkos::subview(neighbor_scratch, team_member.team_rank(), Kokkos::ALL());
//      uint32_t n_neighbors = 0;
//      n_neighbors = EssensKokkos::get_num_neighbors(graph, node, 4, neighbor_subview);
//      num_neighbors(node) = n_neighbors;
//      for(int i=1; i<5; i++) {
//        num_combinations(node) += get_num_combinations(n_neighbors, i);
//      }
//    });
//  });
//printf("Calculated # of combinations\n");
//  for(int i=0; i<num_combinations.extent(0); i++) {
//    printf("%ld ", num_combinations(i));
//  }
//  printf("\n");
//  uint64_t total_combinations = 0;
//  Kokkos::parallel_reduce("Number of combinations", graph.numRows(), 
//  KOKKOS_LAMBDA(const int i, uint64_t& update) {
//    update += num_combinations(i);
//  }, total_combinations);
  uint64_t num_intervals = total_combinations/k_interval;
  if(num_intervals*k_interval < total_combinations) {
    num_intervals++;
  }
  printf("Computed # of combinations: (%lu) split into %d groups\n", total_combinations, num_intervals);
#ifdef DEBUG
  printf("Computed # of combinations: (%lu) split into %d groups\n", total_combinations, num_intervals);
#endif
  Kokkos::View<uint32_t*> starts("Start indices", num_intervals+1);
  starts(0) = 0;
  Kokkos::parallel_for("Find starts", Kokkos::RangePolicy<>(0, starts.extent(0)-1), KOKKOS_LAMBDA(const int i) {
    uint64_t threshold = i*k_interval;
    uint64_t counter = 0;
    for(uint32_t j=0; j<num_combinations.extent(0); j++) {
      counter += num_combinations(j);
      if(counter > threshold) {
        starts(i) = j;
        break;
      }
    }
  });
  starts(starts.extent(0)-1) = graph.numRows();
#ifdef DEBUG
  printf("Rank %d: ", rankn);
  for(uint32_t i=0; i<starts.extent(0); i++) {
    printf("%ld ", starts(i));
  }
  printf("\n");
#endif

  if(comm_size == 1)
  {
    Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Dynamic>> policy(1, NUM_THREADS);
    using member_type = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Dynamic>>::member_type;

#ifdef DEBUG
    cout << "Team size: " << policy.team_size() << endl;
#endif
    Kokkos::View<int**> all_neighbors("Neighbor scratch", policy.team_size(), graph.numRows());
//    Kokkos::Experimental::ScatterView<uint32_t**> metrics_sa(graph_GDV);
    Kokkos::Experimental::ScatterView<uint32_t**, LAYOUT> metrics_sa(graph_GDV);

    Kokkos::View<int**> indices("Index", policy.team_size(), 5);
    Kokkos::View<int**> combination_view("combination", policy.team_size(), 5);
    Kokkos::View<int**> sgraph_distance_signature("dist sig", policy.team_size(), orbits.distance.extent(0));
    Kokkos::View<int**> sgraph_degree_signature("Degree signature", policy.team_size(), 5);
    Kokkos::View<int*[5][5]> induced_subgraph("Subgraph", policy.team_size());
    Kokkos::View<bool** > visited("BFS visited", policy.team_size(), 5);
    Kokkos::View<int**> queue("BFS queue", policy.team_size(), 5);
    Kokkos::View<int**> distance("BFS distance", policy.team_size(), 5);
#ifdef DEBUG
    cout << "Allocated data\n";
#endif

    int i=0;
#ifdef AUTO_CHECKPOINT
    auto ctx = KokkosResilience::make_context(MPI_COMM_SELF, "/home/ntan1/Src_Fido_Kokkos/fido.json");
    printf("Created context\n");
    const auto filt = KokkosResilience::Filter::NthIterationFilter(1);
    printf("Created filter\n");
    i = KokkosResilience::latest_version(*ctx, graph_name);
    if(i < 0)
      i = 0;
    printf("Got latest counter %d\n", i);
#endif
#ifdef ANALYSIS
    int checkpoint_num = 0;
#ifdef HASH_ANALYSIS
    uint64_t page_size = sysconf(_SC_PAGE_SIZE);
//    uint64_t page_size = 2048;
    Kokkos::View<uint64_t*> old_hashes("Old hashes", graph_GDV.span()*sizeof(uint32_t)/page_size);
    Kokkos::View<uint64_t*> new_hashes("New hashes", graph_GDV.span()*sizeof(uint32_t)/page_size);
    auto& old_data = old_hashes;
    auto& new_data = new_hashes;
#endif
#ifdef SCAN_ANALYSIS
    GDVs old_updates("Old updates", graph_GDV.layout());
    auto& old_data = old_updates;
    auto& new_data = graph_GDV;
#endif
#endif
    for(i; i<starts.extent(0)-1; i++) {

#ifdef DIRTY_PAGE_TRACKING
      printf("Chunk %d\n", i);
      reset_dirty_bit(pid);
#endif

#ifdef AUTO_CHECKPOINT
      KokkosResilience::checkpoint(*ctx, graph_name, i, [=] () mutable {
#endif
      chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
      Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Dynamic>> bundle_policy(1, NUM_THREADS);
      Kokkos::parallel_for("Calcualte GDV bundle", bundle_policy, KOKKOS_LAMBDA(member_type team_member) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, starts(i+1)-starts(i)), [=] (int n_offset) {
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
          kokkos_calculate_GDV(team_member, node, graph, orbits, neighbor_subview, indices_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, metrics_sa);
//          kokkos_calculate_GDV(team_member, node, graph, orbits, neighbor_subview, indices_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, graph_GDV);
//          printf("Done with node: %d, time: %f\n", node, time_span.count());
        });
      });
      Kokkos::Experimental::contribute(graph_GDV, metrics_sa);
      metrics_sa.reset();
#ifdef ANALYSIS
std::string region_log("region-data-graph-");
region_log = region_log + std::string(graph_name) + std::string(".log");
std::fstream fs(region_log, std::fstream::out|std::fstream::app);
// Hash based change tracking
uint32_t num_changed = 0;
#ifdef HASH_ANALYSIS
generate_hashes(graph_GDV, new_hashes, page_size);
#endif
//fs << "Number of blocks: " << new_hashes.size() << std::endl;
fs << "Number of blocks: " << new_data.size() << std::endl;
fs << "Changed blocks: ";
bool flag = false;
//for(int idx = 0; idx<new_hashes.size(); idx++) {
//  if(old_hashes(idx) != new_hashes(idx))
for(int idx = 0; idx<new_data.size(); idx++) {
  if(old_data.data()[idx] != new_data.data()[idx])
  {
    num_changed += 1;
    if(!flag) {
      fs << "[" << idx << ",";
      flag = true;
    }
  } else {
    if(flag) {
      fs << idx << ") ";
      flag = false;
    }
  }
}
if(flag)
  fs << new_data.size() << ")";
//  fs << new_hashes.size() << ")";
fs << std::endl;
fs << num_changed << "/" << new_data.size() << " blocks changed " << std::endl;
//fs << num_changed << "/" << new_hashes.size() << " blocks changed " << std::endl;
// Find contiguous regions
std::map<int, int> contiguous_regions;
int largest_region = 1;
int region_size_counter = 0;
//for(int idx=0; idx<new_hashes.size(); idx++) {
//  if(old_hashes(idx) != new_hashes(idx))
for(int idx=0; idx<new_data.size(); idx++) {
  if(old_data.data()[idx] != new_data.data()[idx])
  {
    region_size_counter += 1;
  }
  else 
  {
    if(region_size_counter > largest_region)
    {
      largest_region = region_size_counter;
    }
    if(region_size_counter > 0) {
      auto pos = contiguous_regions.find(region_size_counter);
      if(pos != contiguous_regions.end())
      {
        pos->second = pos->second + 1;
      }
      else
      {
        contiguous_regions.insert(std::pair<int,int>(region_size_counter, 1));
      }
    }
    region_size_counter = 0;
  }
}
if(region_size_counter > largest_region)
{
  largest_region = region_size_counter;
}
if(region_size_counter > 0) {
  auto pos = contiguous_regions.find(region_size_counter);
  if(pos != contiguous_regions.end())
  {
    pos->second = pos->second + 1;
  }
  else
  {
    contiguous_regions.insert(std::pair<int,int>(region_size_counter, 1));
  }
}
fs << "Largest region: " << largest_region << std::endl;
fs << "Region map\n";
for(auto itr = contiguous_regions.begin(); itr != contiguous_regions.end(); ++itr)
{
  fs << itr->second << " regions of size " << itr->first << std::endl;
}
#ifdef HASH_ANALYSIS
//Kokkos::deep_copy(old_hashes, new_hashes);
//Kokkos::deep_copy(new_hashes, 0);
Kokkos::deep_copy(old_data, new_data);
Kokkos::deep_copy(new_data, 0);
#endif
#ifdef SCAN_ANALYSIS
Kokkos::deep_copy(old_data, new_data);
#endif
#endif

#ifdef DIRTY_PAGE_TRACKING
      uintptr_t start_addr = reinterpret_cast<uintptr_t>(graph_GDV.data());
      uintptr_t end_addr = reinterpret_cast<uintptr_t>(graph_GDV.data() + graph_GDV.span());
      bool range_dirty = address_range_dirty(pid, start_addr, end_addr);
      if(range_dirty) {
        printf("Detected changes\n");
      }
      size_t num_pages = ((end_addr-start_addr)/page_size);
      uint64_t* page_list = (uint64_t*) malloc(sizeof(uint64_t)*num_pages);
      for(int j=0; j<num_pages; j++) {
        page_list[j] = 0; 
      }
      bool dirty = get_dirty_pages(pid, start_addr, end_addr, page_list);
      printf("List of dirty pages: ");
      for(int j=0; j<num_pages; j++) {
        if(page_list[j] > 0) {
          printf("%llu, ", page_list[j]);
        }
      }
      printf("\n");
#endif

      chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
      chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2-t1);
//#ifdef DEBUG
      printf("Done with chunk: %d, time: %f\n", i, time_span.count());
//#endif

#ifdef AUTO_CHECKPOINT
}, filt);
#endif
    }
  }
  else
  {
#ifdef ANALYSIS
    int checkpoint_num = 0;
//    uint64_t page_size = sysconf(_SC_PAGE_SIZE);
    uint64_t page_size = 4096;
    Kokkos::View<uint64_t*> old_hashes("Old hashes", graph_GDV.span()*sizeof(uint32_t)/page_size);
    Kokkos::View<uint64_t*> new_hashes("New hashes", graph_GDV.span()*sizeof(uint32_t)/page_size);
#endif
    Kokkos::deep_copy(graph_GDV, 0);
    uint32_t intervals_per_rank = num_intervals/comm_size;
    uint32_t remaining_intervals = num_intervals%comm_size;
    if(rankn < remaining_intervals)
      intervals_per_rank += 1;
    uint32_t start_index = 0;
    if(rankn < remaining_intervals)
      start_index = intervals_per_rank*rankn;
    if(rankn >= remaining_intervals)
      start_index = (intervals_per_rank+1)*remaining_intervals+intervals_per_rank*(rankn-remaining_intervals);
    if(start_index+intervals_per_rank > num_intervals)
      intervals_per_rank = num_intervals - start_index;
printf("Rank %d start index: %d, intervals per rank: %d\n", rankn, start_index, intervals_per_rank);
//    if(rankn == comm_size-1) {
//      intervals_per_rank += num_intervals - (comm_size*intervals_per_rank);
//    }
#ifdef DEBUG
printf("Rank %d start index: %d, intervals per rank: %d\n", rankn, start_index, intervals_per_rank);
#endif
    int offset = 0;
    Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(1, NUM_THREADS);
    using member_type = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;

    Kokkos::Experimental::ScatterView<uint32_t**, LAYOUT> metrics_sa(graph_GDV);
//    Kokkos::Experimental::ScatterView<GDVs::data_type, 
//                                      GDVs::array_layout, 
//                                      GDVs::execution_space, 
//                                      Kokkos::Experimental::ScatterSum, 
//                                      Kokkos::Experimental::ScatterNonDuplicated, 
//                                      Kokkos::Experimental::ScatterAtomic> metrics_sa(graph_GDV);

    Kokkos::View<int**> all_neighbors("Neighbor scratch", policy.team_size(), graph.numRows());
    Kokkos::View<int**> indices("Index", policy.team_size(), 5);
    Kokkos::View<int**> combination_view("combination", policy.team_size(), 5);
    Kokkos::View<int**> sgraph_distance_signature("dist sig", policy.team_size(), orbits.distance.extent(0));
    Kokkos::View<int**> sgraph_degree_signature("Degree signature", policy.team_size(), 5);
    Kokkos::View<int*[5][5]> induced_subgraph("Subgraph", policy.team_size());
    Kokkos::View<bool** > visited("BFS visited", policy.team_size(), 5);
    Kokkos::View<int**> queue("BFS queue", policy.team_size(), 5);
    Kokkos::View<int**> distance("BFS distance", policy.team_size(), 5);
#ifdef DEBUG
printf("Rank %d allocated memory\n", rankn);
#endif

#ifdef AUTO_CHECKPOINT
    auto ctx = KokkosResilience::make_context(MPI_COMM_SELF, "/home/ntan1/Src_Fido_Kokkos/fido.json");
    printf("Created context\n");
    const auto filt = KokkosResilience::Filter::NthIterationFilter(1);
    printf("Created filter\n");
    offset = KokkosResilience::latest_version(*ctx, label);
    if(offset < 0)
      offset = 0;
    printf("Got latest counter %d\n", offset);
#endif



    while(offset < intervals_per_rank) {
//printf("Offset %u\n", offset);
//printf("Rank %d, offset %d\n", rankn, offset);
      chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
      uint32_t chunk_idx = start_index+offset;

#ifdef AUTO_CHECKPOINT
KokkosResilience::checkpoint(*ctx, label, offset, [=] () mutable {
#endif
#ifdef DIRTY_PAGE_TRACKING
      reset_dirty_bit(pid);
#endif
      bool multi_node = (chunk_idx+1 == starts.extent(0) || starts(chunk_idx+1) != starts(chunk_idx)) && 
                        (chunk_idx==0 || starts(chunk_idx-1) != starts(chunk_idx));
      if(multi_node) {
        int start_node = starts(chunk_idx);
        int end_node;
        if(chunk_idx+1 == starts.extent(0)) {
          end_node = graph.numRows();
        } else {
          end_node = starts(chunk_idx+1);
        }
        Kokkos::parallel_for("Calculate GDV", policy, KOKKOS_LAMBDA(member_type team_member) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, end_node-start_node), [=] (int n) {
            int node = start_node + n;
            auto neighbor_subview = Kokkos::subview(all_neighbors, team_member.team_rank(), Kokkos::ALL());
            auto indices_subview = Kokkos::subview(indices, team_member.team_rank(), Kokkos::ALL());
            auto combination_subview = Kokkos::subview(combination_view, team_member.team_rank(), Kokkos::ALL());
            auto sgraph_dist_subview = Kokkos::subview(sgraph_distance_signature, 
                                                        team_member.team_rank(), 
                                                        Kokkos::ALL());
            auto sgraph_deg_subview = Kokkos::subview(sgraph_degree_signature, 
                                                      team_member.team_rank(), 
                                                      Kokkos::ALL());
            auto subgraph_subview = Kokkos::subview(induced_subgraph, 
                                                    team_member.team_rank(), 
                                                    Kokkos::ALL(), 
                                                    Kokkos::ALL());
            auto visited_subview = Kokkos::subview(visited, team_member.team_rank(), Kokkos::ALL());
            auto queue_subview = Kokkos::subview(queue, team_member.team_rank(), Kokkos::ALL());
            auto distance_subview = Kokkos::subview(distance, team_member.team_rank(), Kokkos::ALL());
            kokkos_calculate_GDV(team_member, node, graph, orbits, neighbor_subview, indices_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, metrics_sa);
//            kokkos_calculate_GDV(team_member, node, graph, orbits, neighbor_subview, indices_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, graph_GDV);
          });
        });
#ifdef DEBUG
printf("Rank %d done with chunk %d\n", rankn, chunk_idx);
#endif
      } else {
        auto neighbor_scratch = Kokkos::subview(all_neighbors, 0, Kokkos::ALL());
        int n_neighbors = EssensKokkos::find_neighbours(starts(chunk_idx), graph, 4, neighbor_scratch);

        size_t start_combination;
        size_t end_combination;
        size_t start_comb_subgraph[5] = {0,0,0,0,0};
        size_t end_comb_subgraph[5] = {0,0,0,0,0};
        start_comb_subgraph[0] = 0;
        end_comb_subgraph[0] = 0;
        int start_chunk = 0;
        for(int j=0; j<chunk_idx; j++) {
          if(starts(j) == starts(chunk_idx))
            start_chunk++;
        }
        bool first_chunk=false, middle_chunk=false, last_chunk=false;
        first_chunk = ((chunk_idx == 0) || 
                            starts(chunk_idx-1) != starts(chunk_idx)) && 
                            starts(chunk_idx+1) == starts(chunk_idx);
        if(!first_chunk) {
          last_chunk = starts(chunk_idx-1) == starts(chunk_idx) && 
                            starts(chunk_idx+1) != starts(chunk_idx);
          middle_chunk = starts(chunk_idx-1) == starts(chunk_idx) && 
                              starts(chunk_idx+1) == starts(chunk_idx);
        }
        if(last_chunk) {
          // Last chunk
          start_combination = start_chunk*k_interval;
          end_combination = num_combinations(starts(chunk_idx));
          int64_t counter = end_combination;
          for(int j=4; j>0; j--) {
            int64_t n_comb = get_num_combinations(n_neighbors, j);
            end_comb_subgraph[j] = n_comb;
            counter -= n_comb;
            if(counter > start_combination) {
              start_comb_subgraph[j] = 0;
            } else {
              start_comb_subgraph[j] = start_combination-counter;
              break;
            }
          }
        } else if(first_chunk) {
          // First chunk
          start_combination = 0;
          end_combination = k_interval;
          size_t counter = k_interval;
          for(int j=1; j<5; j++) {
            int64_t n_comb = get_num_combinations(n_neighbors, j);
            if(counter > n_comb) {
              end_comb_subgraph[j] = n_comb;
              counter -= n_comb;
            } else {
              end_comb_subgraph[j] = counter;
              break;
            }
          }
        } else if(middle_chunk) {
          // Middle chunk
          start_combination = start_chunk*k_interval;
          end_combination = start_combination+k_interval;
          size_t counter = 0;
          for(int j=1; j<5; j++) {
            int64_t n_comb = get_num_combinations(n_neighbors, j);
            if(start_combination > counter && start_combination < counter+n_comb) {
              start_comb_subgraph[j] = start_combination-counter;
            }
            if(end_combination > counter+n_comb && counter+n_comb >= start_combination) {
              end_comb_subgraph[j] = n_comb;
            }
            if(end_combination > counter && end_combination < counter+n_comb) {
              end_comb_subgraph[j] = end_combination-counter;
              break;
            }
            counter += n_comb;
          }
        }
        
        int node = starts(chunk_idx);
        Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> team_policy(1, NUM_THREADS);
        Kokkos::View<int**> scratch_block("Indices scratch", team_policy.team_size(), graph.numRows());
        for(int node_count = 1; node_count < 5; node_count++) {
          int64_t s_comb = start_comb_subgraph[node_count];
          int64_t e_comb = end_comb_subgraph[node_count];
          if(e_comb-s_comb > 0) {
            Kokkos::parallel_for("Calculate GDV", team_policy, KOKKOS_LAMBDA(member_type team_member) {
              Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, e_comb-s_comb), [=] (int64_t idx) {
                int64_t combination_num = idx+s_comb;
                auto indices_subview = Kokkos::subview(indices, team_member.team_rank(), Kokkos::ALL());
                auto combination_subview = Kokkos::subview(combination_view, team_member.team_rank(), Kokkos::ALL());
                auto sgraph_dist_subview = Kokkos::subview(sgraph_distance_signature, 
                                                            team_member.team_rank(), 
                                                            Kokkos::ALL());
                auto sgraph_deg_subview = Kokkos::subview(sgraph_degree_signature, 
                                                          team_member.team_rank(), 
                                                          Kokkos::ALL());
                auto subgraph_subview = Kokkos::subview(induced_subgraph, 
                                                        team_member.team_rank(), 
                                                        Kokkos::ALL(), 
                                                        Kokkos::ALL());
                auto visited_subview = Kokkos::subview(visited, team_member.team_rank(), Kokkos::ALL());
                auto queue_subview = Kokkos::subview(queue, team_member.team_rank(), Kokkos::ALL());
                auto distance_subview = Kokkos::subview(distance, team_member.team_rank(), Kokkos::ALL());
                auto scratch_view = Kokkos::subview(scratch_block, team_member.team_rank(), Kokkos::ALL());
                combination_from_position(scratch_view, combination_num, n_neighbors, node_count);
                for(int j=0; j<node_count; j++) {
//                    combination_subview(j) = neighbor_subview(scratch_view(j));
                    combination_subview(j) = neighbor_scratch(scratch_view(j));
                }
                combination_subview(node_count) = node;
//                kokkos_calculate_GDV(team_member, node, node_count, graph, orbits, neighbor_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, metrics_sa);
                kokkos_calculate_GDV(team_member, node, node_count, graph, orbits, neighbor_scratch, num_neighbors(starts(chunk_idx)), combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, metrics_sa);
              });
            });
          }
        }
#ifdef DEBUG
printf("Rank %d done with chunk %d\n", rankn, chunk_idx);
#endif
      }
      Kokkos::Experimental::contribute(graph_GDV, metrics_sa);
      metrics_sa.reset();
      Kokkos::fence();
#ifdef ANALYSIS
std::string region_log("Rank-");
region_log = region_log + std::to_string(rankn) + std::string("-region-data-graph-");
region_log = region_log + std::string(graph_name) + std::string(".log");
std::fstream fs(region_log, std::fstream::out|std::fstream::app);
// Hash based change tracking
uint32_t num_changed = 0;
generate_hashes(graph_GDV, new_hashes, page_size);
fs << "Number of blocks: " << new_hashes.size() << std::endl;
fs << "Changed blocks: ";
bool flag = false;
for(int idx = 0; idx<new_hashes.size(); idx++) {
  if(old_hashes(idx) != new_hashes(idx))
  {
    num_changed += 1;
    if(!flag) {
      fs << "[" << idx << ",";
      flag = true;
    }
  } else {
    if(flag) {
      fs << idx << ") ";
      flag = false;
    }
  }
}
if(flag)
  fs << new_hashes.size() << ")";
//std::cout << std::endl;
//std::cout << num_changed << "/" << new_hashes.size() << " blocks changed " << std::endl;
fs << std::endl;
fs << num_changed << "/" << new_hashes.size() << " blocks changed " << std::endl;
// Find contiguous regions
std::map<int, int> contiguous_regions;
int largest_region = 1;
int region_size_counter = 0;
for(int idx=0; idx<new_hashes.size(); idx++) {
  if(old_hashes(idx) != new_hashes(idx))
  {
    region_size_counter += 1;
  }
  else 
  {
    if(region_size_counter > largest_region)
    {
      largest_region = region_size_counter;
    }
    if(region_size_counter > 0) {
      auto pos = contiguous_regions.find(region_size_counter);
      if(pos != contiguous_regions.end())
      {
        pos->second = pos->second + 1;
      }
      else
      {
        contiguous_regions.insert(std::pair<int,int>(region_size_counter, 1));
      }
    }
    region_size_counter = 0;
  }
}
if(region_size_counter > largest_region)
{
  largest_region = region_size_counter;
}
if(region_size_counter > 0) {
  auto pos = contiguous_regions.find(region_size_counter);
  if(pos != contiguous_regions.end())
  {
    pos->second = pos->second + 1;
  }
  else
  {
    contiguous_regions.insert(std::pair<int,int>(region_size_counter, 1));
  }
}
//std::cout << "Largest region: " << largest_region << std::endl;
//std::cout << "Region map\n";
fs << "Largest region: " << largest_region << std::endl;
fs << "Region map\n";
for(auto itr = contiguous_regions.begin(); itr != contiguous_regions.end(); ++itr)
{
//  std::cout << itr->second << " regions of size " << itr->first << std::endl;
  fs << itr->second << " regions of size " << itr->first << std::endl;
}
Kokkos::deep_copy(old_hashes, new_hashes);
Kokkos::deep_copy(new_hashes, 0);
#endif
#ifdef DIRTY_PAGE_TRACKING
      uintptr_t start_addr = reinterpret_cast<uintptr_t>(graph_GDV.data());
      uintptr_t end_addr = reinterpret_cast<uintptr_t>(graph_GDV.data() + graph_GDV.span());
      bool range_dirty = address_range_dirty(pid, start_addr, end_addr);
      if(range_dirty) {
        printf("Detected changes\n");
      }
      size_t num_pages = ((end_addr-start_addr)/page_size);
      uint64_t* page_list = (uint64_t*) malloc(sizeof(uint64_t)*num_pages);
      for(int j=0; j<num_pages; j++) {
        page_list[j] = 0; 
      }
      bool dirty = get_dirty_pages(pid, start_addr, end_addr, page_list);
      printf("List of dirty pages: ");
      for(int j=0; j<num_pages; j++) {
        if(page_list[j] > 0) {
          printf("%llu, ", page_list[j]);
        }
      }
      printf("\n");
#endif
#ifdef AUTO_CHECKPOINT
}, filt);
#endif
      chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
      chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2-t1);
      printf("Done with chunk: %d, time: %f\n", chunk_idx, time_span.count());
      offset++;
    }
    printf("Rank %d completed calculation\n", rankn);
    if(rankn == 0) {
      MPI_Reduce(MPI_IN_PLACE, graph_GDV.data(), graph_GDV.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
      MPI_Reduce(graph_GDV.data(), graph_GDV.data(), graph_GDV.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //vec_calc_post_gather = MPI_Wtime() - vec_calc_start + vec_calc_post_gather;
  vec_calc_communication_time[graph_counter - 1] = process_ends_communication - vec_calc_start;
//  cout << "Communication time on process " << rankn << " for graph " << graph_counter << " = " << vec_calc_communication_time[graph_counter - 1] << endl;
  //vec_calc_computation_time = vec_calc_computation_end - vec_calc_computation_start + vec_calc_computation_time;

  #ifdef DEBUG
    cout << "Finished GDV Vector Calc on Rank: " << rankn << endl;
  #endif

  printf("Rank %d finished GDV vector calc\n", rankn);
}

template<class NeighborView, class IntView, class GraphView, class BoolView, class Scatter>
KOKKOS_INLINE_FUNCTION void 
kokkos_calculate_GDV(Kokkos::TeamPolicy<>::member_type team_member,
                      int node, 
                      int node_count,
                      const matrix_type& graph, 
                      const Orbits& orbits, 
                      const NeighborView& neighbors,
                      const int num_neighbors,
//                      const IntView& indices,
                      const IntView& combination_view,
                      IntView& sgraph_distance_signature,
                      IntView& sgraph_degree_signature,
                      GraphView& induced_subgraph,
                      BoolView& visited,
                      IntView& queue,
                      IntView& distance,
//                      Kokkos::Experimental::ScatterView<int**, Kokkos::Experimental::ScatterNonDuplicated> gdvMetrics_sa
//                      Kokkos::Experimental::ScatterView<int**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace, Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated, Kokkos::Experimental::ScatterAtomic> gdvMetrics_sa
                      Scatter gdvMetrics_sa
                    )
{
  auto gdvMetrics = gdvMetrics_sa.access();
//  int num_neighbors = neighbors.extent(0);
//  int num_neighbors = EssensKokkos::find_neighbours(node, graph, 4, neighbor_buff);
//  auto neighbors = Kokkos::subview(neighbor_buff, std::pair<int,int>(0, num_neighbors));
//  printf("Allocated and found neighbors of %d: %d total neighbors\n", node, neighbors.size());
//  for (int node_count = 1; node_count < 5; node_count++)
//  {
//    CombinationGenerator generator(num_neighbors, node_count, indices);
//    while(!generator.done) {
      auto combination = Kokkos::subview(combination_view, std::pair<int,int>(0,node_count+1));
//      kokkos_get_combination(indices, node_count, neighbors, combination);
//      for(int i=0; i<node_count; i++) {
//          combination(i) = neighbors(indices(i));
//      }
//      combination(node_count) = node;
      auto subgraph_degree_signature = Kokkos::subview(sgraph_degree_signature, 
                                                        std::pair<int,int>(0, node_count+1));
      auto subgraph_distance_signature = Kokkos::subview(sgraph_distance_signature, 
                                                        std::pair<int,int>(0, orbits.distance.extent(1)));
      auto induced_sgraph = Kokkos::subview(induced_subgraph, 
                                            std::pair<int,int>(0,node_count+1), 
                                            std::pair<int,int>(0,node_count+1));
      EssensKokkos::kokkos_induced_subgraph(graph, combination, induced_sgraph);
      bool is_connected = EssensKokkos::is_connected(induced_sgraph, visited, queue);
      if(is_connected)
      {
        EssensKokkos::calc_degree_signature(induced_sgraph, subgraph_degree_signature);
        for(int idx=0; idx<node_count+1; idx++)
        {
          int v = idx;
          EssensKokkos::calc_distance_signature(v, induced_sgraph, subgraph_distance_signature, 
                                                visited, queue, distance);
          for(int i=orbits.start_indices(node_count+1); i<orbits.start_indices(node_count+2); i++) {
            auto orbit_deg_sig = Kokkos::subview(orbits.degree, i, Kokkos::ALL);
            bool match = EssensKokkos::compare_signatures(subgraph_degree_signature, orbit_deg_sig);
            auto orbit_dis_sig = Kokkos::subview(orbits.distance, i, Kokkos::ALL);
            match = match && EssensKokkos::compare_signatures(subgraph_distance_signature, orbit_dis_sig);
            if(match) {
              gdvMetrics(combination(v),i) += 1;
            }
          }
        }
      }
//      generator.kokkos_next(indices);
//    }
//    printf("Number of combinations: %d for %d neighbors and %d node subgraphs\n", generator.get_num_comb(), num_neighbors, node_count+1);
//  }
//  counter(team_member.team_rank()) = iter_counter;
//  printf("Thread %d at iteration %d\n", team_member.team_rank(), iter_counter);
}


//template<class NeighborView, class IntView, class GraphView, class BoolView, class CounterView>
template<class NeighborView, class IntView, class GraphView, class BoolView, class Scatter>
KOKKOS_INLINE_FUNCTION void 
kokkos_calculate_GDV(Kokkos::TeamPolicy<>::member_type team_member,
                      int node, 
                      const matrix_type& graph, 
                      const Orbits& orbits, 
                      NeighborView& neighbor_buff,
                      IntView& indices,
                      IntView& combination_view,
                      IntView& sgraph_distance_signature,
                      IntView& sgraph_degree_signature,
                      GraphView& induced_subgraph,
                      BoolView& visited,
                      IntView& queue,
                      IntView& distance,
//                      CounterView& counter,
//                      Kokkos::Experimental::ScatterView<int**, Kokkos::Experimental::ScatterNonDuplicated> gdvMetrics_sa
                      Scatter gdvMetrics_sa
                    )
{
  auto gdvMetrics = gdvMetrics_sa.access();
//  int combination_count = 0;
//  int k_interval = 1000000;
  int num_neighbors = EssensKokkos::find_neighbours(node, graph, 4, neighbor_buff);
  auto neighbors = Kokkos::subview(neighbor_buff, std::pair<int,int>(0, num_neighbors));
//  printf("Allocated and found neighbors of %d: %d total neighbors\n", node, neighbors.size());
//  int iter_counter = counter(team_member.team_rank());
  for (int node_count = 1; node_count < 5; node_count++)
  {
    CombinationGenerator generator(num_neighbors, node_count, indices);
    while(!generator.done) {
      auto combination = Kokkos::subview(combination_view, std::pair<int,int>(0,node_count+1));
      generator.kokkos_get_combination(indices, num_neighbors, neighbors, combination);
      auto subgraph_degree_signature = Kokkos::subview(sgraph_degree_signature, 
                                                        std::pair<int,int>(0, node_count+1));
      auto subgraph_distance_signature = Kokkos::subview(sgraph_distance_signature, 
                                                        std::pair<int,int>(0, orbits.distance.extent(1)));
      combination(node_count) = node;
      auto induced_sgraph = Kokkos::subview(induced_subgraph, 
                                            std::pair<int,int>(0,node_count+1), 
                                            std::pair<int,int>(0,node_count+1));
      EssensKokkos::kokkos_induced_subgraph(graph, combination, induced_sgraph);
      bool is_connected = EssensKokkos::is_connected(induced_sgraph, visited, queue);
      if(is_connected)
      {
        EssensKokkos::calc_degree_signature(induced_sgraph, subgraph_degree_signature);
        for(int idx=0; idx<node_count+1; idx++)
        {
          int v = idx;
          EssensKokkos::calc_distance_signature(v, induced_sgraph, subgraph_distance_signature, 
                                                visited, queue, distance);
          for(int i=orbits.start_indices(node_count+1); i<orbits.start_indices(node_count+2); i++) {
            auto orbit_deg_sig = Kokkos::subview(orbits.degree, i, Kokkos::ALL);
            bool match = EssensKokkos::compare_signatures(subgraph_degree_signature, orbit_deg_sig);
            auto orbit_dis_sig = Kokkos::subview(orbits.distance, i, Kokkos::ALL);
            match = match && EssensKokkos::compare_signatures(subgraph_distance_signature, orbit_dis_sig);
            if(match) {
//              Kokkos::atomic_increment(&gdvMetrics_sa(combination(v),i));
              gdvMetrics(combination(v),i) += 1;
//if(node == 1) {
//if(comb_counter < 100)
//  n_comb++;
//}
  
            }
          }
        }
      }
      generator.kokkos_next(indices);
//if(node == 1) {
//if(comb_counter < 100) {
//  printf("Number of updates for iter %d: %d\n", comb_counter, n_comb);
//  n_comb = 0;
//  comb_counter++;
//}
//}

//      combination_count++;
//      iter_counter += 1;
    }
//    printf("Number of combinations: %d for %d neighbors and %d node subgraphs\n", generator.get_num_comb(), num_neighbors, node_count+1);
  }
//  counter(team_member.team_rank()) = iter_counter;
//  printf("Thread %d at iteration %d\n", team_member.team_rank(), iter_counter);
}

void Calculate_GDV(int node,A_Network Graph,vector<OrbitMetric> &orbits, vector<GDVMetric> &gdvMetrics)
{
  GDV_functions gdvf;
  vector<int> neighbours;
  gdvf.find_neighbours(node,Graph,4,&neighbours);
  int set[neighbours.size()]; 
  std::copy( neighbours.begin(), neighbours.end(), set );
//  int numElements = sizeof(set)/sizeof(set[0]);
  for (int node_count = 1; node_count < 5; node_count++)
  {
//    vector<vector<int>> combinationsList;
//    gdvf.find_combinations(set, numElements,node_count,&combinationsList);
    CombinationGeneratorVector generator(neighbours.size(), node_count);
    vector<int> combination(node_count, 0);
//    for (vector<int> &combination : combinationsList)
//    {
    while(!generator.done) {
      generator.get_combination(neighbours, combination);
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
        sort(subgraph_degree_signature.begin(), subgraph_degree_signature.end());
        vector<OrbitMetric> filter_orbits;
        gdvf.orbit_filter(orbits, node_count+1, filter_orbits);
        for(size_t idx=0; idx<induced_sgraph.size(); idx++)
        {
          int v = induced_sgraph[idx].Row;
          gdvf.distance_signature(v,induced_sgraph,subgraph_distance_signature);
          for(OrbitMetric orbit: filter_orbits)
          {
            if( orbit.orbitDistance == subgraph_distance_signature && 
                orbit.orbitDegree == subgraph_degree_signature)
            {
              gdvMetrics[v].GDV[orbit.orbitNumber] += 1;
//              break;
            }
          }
        }
      }
      generator.next();
    }
  }
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

    sort(vector_line[1].begin(), vector_line[1].end());
    OrbitMetric orbMetric(vector_line[0][0],vector_line[1],vector_line[2]);
    output->push_back(orbMetric);
  }
}

//This method takes the file and converts it into orbits and saves in output
KOKKOS_INLINE_FUNCTION void kokkos_readin_orbits(ifstream *file, Orbits& orbits )
{
  string line;
  string signature_delimiter;
  string internal_delimiter;
  signature_delimiter = "/";
  internal_delimiter= ",";
  int num_orbits = 0;
  file->clear();
  file->seekg(0);
  while(std::getline(*file,line))
  {
    num_orbits++;
  }
  file->clear();
  file->seekg(0);
  orbits = Orbits(num_orbits, 5);
  orbits.start_indices(0) = 0;
  orbits.start_indices(1) = 0;
  orbits.start_indices(2) = 0;
  size_t orbit_size = 2;
  int orbit_counter = 0;
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
    sort(vector_line[1].begin(), vector_line[1].end());
    if(vector_line[1].size() > orbit_size) {
      orbit_size = vector_line[1].size();
      orbits.start_indices(orbit_size) = orbit_counter;
    }
    for(size_t i=0; i<vector_line[1].size(); i++) {
      orbits.degree(orbit_counter, i) = vector_line[1][i];
    }
    for(size_t i=0; i<vector_line[2].size(); i++) {
      orbits.distance(orbit_counter, i) = vector_line[2][i];
    }
    orbit_counter++;
  }
  orbits.start_indices(6) = num_orbits;
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

bool sort_by_leading_edge(vector<int>& a, vector<int>& b) {
  return a[0] < b[0];
}

bool sort_by_dest(vector<int>& a, vector<int>& b) {
  return a[1] < b[1];
}

KOKKOS_INLINE_FUNCTION void readin_graph(ifstream* file, matrix_type& graph) 
{
  string line;;

  int lastrow = -1;
  vector<int> rows, cols;
  vector<float> vals;
  vector<vector<int>> edge_list;
  while(std::getline(*file,line))
  {
    int u, v;
    float w;
    sscanf(line.c_str(), "%d %d %f", &u, &v, &w);
    if(u > lastrow) 
      lastrow = u;
    if(v > lastrow) 
      lastrow = v;
    vector<int> edge1;
    edge1.push_back(u);
    edge1.push_back(v);
    if(w == 0.0 || w == 0) {
        edge1.push_back(0);
    } else {
        edge1.push_back(1);
    }
    edge_list.push_back(edge1);
    if(u != v) {
      vector<int> edge2;
      edge2.push_back(v);
      edge2.push_back(u);
      if(w == 0.0 || w == 0) {
            edge2.push_back(0);
      } else {
            edge2.push_back(1);
      }
      edge_list.push_back(edge2);
    }
  }
  sort(edge_list.begin(), edge_list.end(), sort_by_leading_edge);

  vector<int> rowmap(lastrow+2, 0);
  for(size_t i=0; i<edge_list.size(); i++) {
    rows.push_back(edge_list[i][0]);
    cols.push_back(edge_list[i][1]);
    vals.push_back(edge_list[i][2]);
  }
  
  for(size_t i=0; i<rows.size(); i++) {
    if(rows[i] != cols[i])
    {
      rowmap[rows[i]+1]++;
    }
    else
    {
      rowmap[rows[i]+1]++;
    }
  }
  for(size_t i=1; i<rowmap.size(); i++) {
    rowmap[i] += rowmap[i-1];
  }
  for(size_t i=0; i<rowmap.size()-1; i++) {
    sort(cols.begin()+rowmap[i], cols.begin()+rowmap[i+1]);
  }
  
  graph = matrix_type("Graph", lastrow+1, lastrow+1, vals.size(), 
                      vals.data(), rowmap.data(), cols.data());
  return;
}

