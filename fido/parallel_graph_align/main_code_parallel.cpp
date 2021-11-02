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
#define NUM_THREADS 1
#define MAX_SUBGRAPH 5
//#define DEBUG
//#define INSTRUMENT
//#define DIRTY_PAGE_TRACKING
//#define RESILIENCE
#define AUTO_CHECKPOINT
//#define STANDARD
//#define INCREMENTAL
//#define HASH_DETECT
//#define OUTPUT_MATRIX
#define OUTPUT_GDV
//#define OUTPUT_MATRIX

#ifdef AUTO_CHECKPOINT
#include <resilience/Resilience.hpp>
#include <resilience/CheckpointFilter.hpp>
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
  
    if (myfile.is_open()) {
      myfile << month << "/" << day << "/" << year << "\t" << numtasks << "\t" << NUM_THREADS << "\t" << graph_name1 << "\t" << graph_name2 << "\t" << graphx.numRows() << "\t\t" << total_time_taken << " \n";
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
  kokkos_GDV_vector_calculation(graph1, graph1_GDV, orbits, "graph1", graph_counter); 
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

  graph_counter = 2;
  kokkos_GDV_vector_calculation(graph2, graph2_GDV, orbits, "graph2", graph_counter); 
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

//void 
//kokkos_GDV_vector_calculation(const matrix_type& graph, 
//                              GDVs& graph_GDV, 
//                              const Orbits& orbits, 
//                              const char* graph_name, 
//                              int graph_counter) {
//  // Set up parallelization               
//  int comm_size, rankn;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rankn);
//  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
//  int tag = 11;
//  int graph_size = graph.numRows();
//cout << "Starting GDV vector calculation\n";
//cout << "# of processes : " << comm_size << endl;
////  int graph_size = graph.size();
////  graph_GDV->clear();
////
////  for(int i=0; i<graph.size(); i++) {
////    vector<int> gdv_vec(GDV_LENGTH, 0);
////    graph_GDV->push_back(GDVMetric(graph[i].Row, gdv_vec));
////  }
//  //graph_counter += 1;
//  //cout << "Graph counter on rank " << rankn << " = " << graph_counter << endl;
//
//  double process_ends_communication;
//  double vec_calc_computation_start;
//  double vec_calc_computation_end;
//  double vec_calc_start = MPI_Wtime();
//
//  Kokkos::View<int*> num_combinations("Number of combinations", graph.numRows());
//  int k_interval = 1000000;
//
//  Kokkos::View<int*> neighbor_scratch("Neighbors", graph.numRows());
//  Kokkos::parallel_for(graph.numRows(), KOKKOS_LAMBDA(const int node) {
//    int num_neighbors = 0;
//    num_neighbors = EssensKokkos::get_num_neighbors(graph, node, 4, neighbor_scratch);
//    for(int i=1; i<5; i++) {
//      num_combinations(node) += get_num_combinations(num_neighbors, i);
//    }
//  });
//  printf("Computed # of combinations\n");
//  Kokkos::parallel_scan(num_combinations.extent(0), KOKKOS_LAMBDA(const int i, int& update, const bool final) {
//    const int val_i = num_combinations(i);
//    if(final) {
//      num_combinations(i) = update;
//    }
//    update += val_i;
//  });
//  Kokkos::View<int*> starts("Start indices", (num_combinations(graph.numRows()-1)/k_interval)+1);
//  Kokkos::View<int*> ends("Start indices", (num_combinations(graph.numRows()-1)/k_interval)+1);
//  starts(0) = 0;
//  int threshold = k_interval;
//  int start_index = 1;
//  for(int i=0; i<num_combinations.extent(0); i++) {
//    if(num_combinations(i) > threshold) {
//      starts(start_index++) = i;
//      threshold += k_interval;
//    }
//  }
//  for(int i=0; i<starts.extent(0)-1; i++) {
//    ends(i) = starts(i+1);
//  }
//  ends(ends.extent(0)-1) = graph.numRows();
//  for(int i=0; i<starts.extent(0); i++) {
//    printf("%d ", starts(i));
//  }
//  printf("\n");
//  for(int i=0; i<starts.extent(0); i++) {
//    printf("%d ", ends(i));
//  }
//  printf("\n");
//  for(int i=0; i<starts.extent(0); i++) {
//    printf("%d ", num_combinations(starts(i)));
//  }
//  printf("\n");
//
//  if(comm_size == 1)
//  {
//#ifdef INSTRUMENT
////printf("Combination | # of updated regions | avg size of region | (# of identical regions, total size)\n");
////printf("Iteration|# of regions|Total size|Avg size of region|%% of entries updated|(updates,empty)|# nonzero|# contiguous regions\n");
//#endif
////    Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>> policy(1, Kokkos::AUTO());
//    Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>> policy(1, NUM_THREADS);
//    using member_type = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>>::member_type;
//
//cout << "Team size: " << policy.team_size() << endl;
//    Kokkos::View<int**> all_neighbors("Neighbor scratch", policy.team_size(), graph.numRows());
//    Kokkos::View<int*> combination_counter("Per thread counter", policy.team_size());
//    Kokkos::deep_copy(combination_counter, 0);
//    Kokkos::Experimental::ScatterView<int**> metrics_sa(graph_GDV);
//
//    Kokkos::View<int**> indices("Index", policy.team_size(), 5);
//    Kokkos::View<int**> combination_view("combination", policy.team_size(), 5);
//    Kokkos::View<int**> sgraph_distance_signature("dist sig", policy.team_size(), orbits.distance.extent(0));
//    Kokkos::View<int**> sgraph_degree_signature("Degree signature", policy.team_size(), 5);
//    Kokkos::View<int*[5][5]> induced_subgraph("Subgraph", policy.team_size());
//    Kokkos::View<bool** > visited("BFS visited", policy.team_size(), 5);
//    Kokkos::View<int**> queue("BFS queue", policy.team_size(), 5);
//    Kokkos::View<int**> distance("BFS distance", policy.team_size(), 5);
//
//    int i=0;
//#ifdef AUTO_CHECKPOINT
//    auto ctx = KokkosResilience::make_context(MPI_COMM_WORLD, "/home/ntan1/Src_Fido_Kokkos/fido.json");
//    printf("Created context\n");
//    const auto filt = KokkosResilience::Filter::NthIterationFilter(1);
//    printf("Created filter\n");
//    i = KokkosResilience::latest_version(*ctx, graph_name);
//    if(i < 0)
//      i = 0;
//    printf("Got latest counter %d\n", i);
//#endif
//    for(i; i<starts.extent(0); i++) {
//#ifdef AUTO_CHECKPOINT
//KokkosResilience::checkpoint(*ctx, graph_name, i, [=] () mutable {
//#endif
//      Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Dynamic>> bundle_policy(1, NUM_THREADS);
//      Kokkos::parallel_for("Calcualte GDV bundle", bundle_policy, KOKKOS_LAMBDA(member_type team_member) {
//        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, ends(i)-starts(i)), [=] (int n_offset) {
//          int node = starts(i) + n_offset;
//          auto neighbor_subview = Kokkos::subview(all_neighbors, team_member.team_rank(), Kokkos::ALL());
//          auto indices_subview = Kokkos::subview(indices, team_member.team_rank(), Kokkos::ALL());
//          auto combination_subview = Kokkos::subview(combination_view, team_member.team_rank(), Kokkos::ALL());
//          auto sgraph_dist_subview = Kokkos::subview(sgraph_distance_signature, team_member.team_rank(), Kokkos::ALL());
//          auto sgraph_deg_subview = Kokkos::subview(sgraph_degree_signature, team_member.team_rank(), Kokkos::ALL());
//          auto subgraph_subview = Kokkos::subview(induced_subgraph, team_member.team_rank(), Kokkos::ALL(), Kokkos::ALL());
//          auto visited_subview = Kokkos::subview(visited, team_member.team_rank(), Kokkos::ALL());
//          auto queue_subview = Kokkos::subview(queue, team_member.team_rank(), Kokkos::ALL());
//          auto distance_subview = Kokkos::subview(distance, team_member.team_rank(), Kokkos::ALL());
//chrono::  steady_clock::time_point t1 = chrono::steady_clock::now();
//          kokkos_calculate_GDV(team_member, node, graph, orbits, neighbor_subview, indices_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, combination_counter, metrics_sa);
//          chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
//          chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2-t1);
//          printf("Done with node: %d, time: %f\n", node, time_span.count());
//        });
//      });
//      Kokkos::Experimental::contribute(graph_GDV, metrics_sa);
//      metrics_sa.reset();
//#ifdef AUTO_CHECKPOINT
//}, filt);
//#endif
//    }
////    Kokkos::parallel_for("Calculate GDV", policy, KOKKOS_LAMBDA(member_type team_member) {
////      Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, graph.numRows()), [=] (int node) {
//////        int node = team_member.league_rank()*team_member.team_size() + team_member.team_rank();
////        auto neighbor_subview = Kokkos::subview(all_neighbors, team_member.team_rank(), Kokkos::ALL());
////        auto indices_subview = Kokkos::subview(indices, team_member.team_rank(), Kokkos::ALL());
////        auto combination_subview = Kokkos::subview(combination_view, team_member.team_rank(), Kokkos::ALL());
////        auto sgraph_dist_subview = Kokkos::subview(sgraph_distance_signature, team_member.team_rank(), Kokkos::ALL());
////        auto sgraph_deg_subview = Kokkos::subview(sgraph_degree_signature, team_member.team_rank(), Kokkos::ALL());
////        auto subgraph_subview = Kokkos::subview(induced_subgraph, team_member.team_rank(), Kokkos::ALL(), Kokkos::ALL());
////        auto visited_subview = Kokkos::subview(visited, team_member.team_rank(), Kokkos::ALL());
////        auto queue_subview = Kokkos::subview(queue, team_member.team_rank(), Kokkos::ALL());
////        auto distance_subview = Kokkos::subview(distance, team_member.team_rank(), Kokkos::ALL());
////chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
////        kokkos_calculate_GDV(team_member, node, graph, orbits, neighbor_subview, indices_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, combination_counter, 
////#ifdef CHECKPOINT_TEST
////                              graph_GDV, chkpt1, chkpt2, 
////#endif
////                              metrics_sa);
////chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
////chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2-t1);
////        printf("Done with node: %d, time: %f\n", node, time_span.count());
////      });
////    });
////    Kokkos::Experimental::contribute(graph_GDV, metrics_sa);
//  }
//  else if (rankn == 0) 
//  {
//    int i;
//    if (graph_size < comm_size) {
//
//      // Send all nodes if comm size is bigger
//      for (i = 0; i < graph_size; i++) {
//        MPI_Send(&i, 1, MPI_INT, i+1, tag, MPI_COMM_WORLD);
//      }
//
//      // Send termination to finish processes
//      int flag;
//      for (i = 1; i < comm_size; i++) {
//        flag = -1;
//        MPI_Send(&flag, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
//      }
//
//    } else { // There are more nodes in graph than there are MPI processes
//
//      // First get each process busy
//      int send_node; // Corresponds to index of node to send
//      int rcv_node;  // Corresponds to name of graph node recieved from worker rank
//      for (i = 1; i < comm_size; i++) {
//	      send_node = (i-1) * CHUNK_SIZE;
//        #ifdef DEBUG
//          cout << "Sending node " << send_node << " to rank " << i << endl;
//        #endif
//        MPI_Send(&send_node, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
//      }
////cout << "Sent initial batch of nodes\n";
//
//      // Start probing and recieving results from other processes
//      int rec_count = 0;
//      send_node += CHUNK_SIZE;
//      //int next_job = comm_size-1;
//      do {
//	
//        // First probe for completed work
//      	int flag;
//      	MPI_Status master_status;
//      	MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &master_status);
//      
//      	if (flag == 1) {
//      	
//      	  // Recieve gdv from finished process
//      	  i = master_status.MPI_SOURCE;
//          Kokkos::View<int[GDV_LENGTH+1]> gdv_array_d("Recv buffer");
//          auto gdv_array_h = Kokkos::create_mirror_view(gdv_array_d);
//      	  MPI_Recv(gdv_array_h.data(), GDV_LENGTH + 1, MPI_INT, i, tag, MPI_COMM_WORLD, &master_status);
//          Kokkos::deep_copy(gdv_array_d, gdv_array_h);
////cout << "Received GDV from " << gdv_array_h(GDV_LENGTH) << endl;
//      	  #ifdef DEBUG
//            cout << "Recieved GDV for node " << rcv_node << " from rank " << i << ": " << endl;
//            cout << gdv_array_h(GDV_LENGTH) << ": ";
//            for (int j = 0; j < GDV_LENGTH; j++) {
//              cout << gdv_array_h(j) << ", ";
//            }
//            cout << endl;
//          #endif
////cout << "Receive buffer for node " << gdv_array_h(GDV_LENGTH) << ": ";
////for(int k=0; k<gdv_array_h.extent(0); k++) {
////  cout << gdv_array_h(k) << " ";
////}
////cout << endl;
//      
//      	  // Organize recieved data into returning array
////      	  rcv_gdv.clear();
////          rcv_gdv.resize(GDV_LENGTH);
////          for (int j = 0; j < GDV_LENGTH; j++) {
////            rcv_gdv[j] = gdv_array[j];
////          }
////          Kokkos::parallel_for("Copy GDV", Kokkos::RangePolicy<>(0, GDV_LENGTH), KOKKOS_LAMBDA(const int j) {
////            rcv_gdv(j) = gdv_array_d[j];
////          });
//          rcv_node = gdv_array_h(GDV_LENGTH);
//          // We are updating each vertex in the nodes neighborhood.
//          // Last GDV update for the node will be node+graph_size 
//          // as a signal for when we're done with this node
//          if(rcv_node >= graph_size) { // Check if last GDV for the node
//            rcv_node -= graph_size;
//            rec_count += CHUNK_SIZE;
////cout << "Received GDVs for node " << rcv_node << endl;
//      	    // Prepare to send next node to finished process.
//      	    if (send_node < graph_size) { // Jobs still exist.  Send next.
//              #ifdef DEBUG
////      	        cout << "Sending node " << graph[send_node].Row << " to rank " << i << endl;
//      	        cout << "Sending node " << send_node << " to rank " << i << endl;
//      	      #endif
//      	      MPI_Send(&send_node, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
//      	      send_node += CHUNK_SIZE;
//      	    } else { // Send termination
//      	      flag = -1;
//      	      MPI_Send(&flag, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
//      	    }
////cout << "Send new batch of nodes to " << i << endl;
//          }
////          for(int j=0; j<GDV_LENGTH; j++) {
////            (*graph_GDV)[rcv_node].GDV[j] += rcv_gdv[j];
////          }
//          for(int j=0; j<GDV_LENGTH; j++) {
//            graph_GDV(rcv_node, j) += gdv_array_d(j);
//          }
////cout << "Updated GDV for " << rcv_node << endl;
////
////cout << "GDV: Root node " << rcv_node << endl;
////for(int i=0; i<graph_GDV.extent(0); i++) {
////  cout << "Node " << i << ": ";
////  for(int j=0; j<graph_GDV.extent(1); j++) {
////    cout << graph_GDV(i,j) << " ";
////  }
////  cout << endl;
////}
////cout << endl;
////      	  free(gdv_array);
//      	}
//      } while (rec_count < graph_size);
//
//      process_ends_communication = MPI_Wtime();
//      //vec_calc_prior_gather = MPI_Wtime() - vec_calc_start + vec_calc_prior_gather;
//
//      // Sort return vector
////      sort(graph_GDV->begin(), graph_GDV->end(), node_id_order());
//      #ifdef DEBUG
//        cout << "Constructed return GDV array" << endl;
//        for (i = 0; i < graph_GDV->size(); i++) {
//      	  cout << graph_GDV->at(i).node << ": ";
//      	  for (int j = 0; j < graph_GDV->at(i).GDV.size(); j++) {
//      	    cout << graph_GDV->at(i).GDV[j] << ", ";
//      	  }
//      	  cout << endl;
//      	}	
//      #endif
//    }
//  }
//  else // Instructions for work processes
//  { 
//    int node_name;
//    MPI_Status worker_status;
//
//    do {
//
////cout << "Starting worker thread\n";
//      MPI_Recv(&node_name, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &worker_status);
////cout << "Got new batch of nodes starting at " << node_name << "\n";
//      #ifdef DEBUG
//        cout << "Recieved node " << node_name << " at rank " << rankn << endl;
//      #endif
//      if (node_name == -1) {
//        process_ends_communication = MPI_Wtime();
//        //vec_calc_prior_gather = MPI_Wtime() - vec_calc_start + vec_calc_prior_gather;
////cout << "Thread " << rankn << " done\n";
//        break;
//      }
//      int end = node_name+CHUNK_SIZE;
////      if(end > graph.size())
////        end = graph.size();
//      int nrows = graph.numRows();
//      if(end > nrows)
//        end = nrows;
////      vector<GDVMetric> metrics;
////      for(int idx=0; idx<graph.numRows(); idx++) {
////        metrics.push_back(GDVMetric(graph[idx].Row, vector<int>(GDV_LENGTH, 0)));
////      }
//#ifdef INSTRUMENT
////printf("Combination | # of updated regions | avg size of region | (# of identical regions, total size)\n");
////printf("Iteration|# of regions|Total size|Avg size of region|%% of entries updated|(updates,empty)|# nonzero|# contiguous regions\n");
//#endif
//
////      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(1, Kokkos::AUTO());
////      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(end-node_name, Kokkos::AUTO());
//      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(1, NUM_THREADS);
////      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy((end-node_name)/4, 4);
////      policy.set_scratch_size(0, Kokkos::PerThread(128));
//      using member_type = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
//
//cout << "Team size: " << policy.team_size() << endl;
//      Kokkos::View<int**> all_neighbors("Neighbor scratch", policy.team_size(), graph.numRows());
//      GDVs metrics("GDVs", graph.numRows(), GDV_LENGTH);
//      Kokkos::Experimental::ScatterView<int**> metrics_sa(metrics);
//
//  Kokkos::View<int*> combination_counter("Per thread counter", policy.team_size());
//  Kokkos::deep_copy(combination_counter, 0);
//  Kokkos::View<int** > indices("Index", policy.team_size(), 5);
//  Kokkos::View<int** > combination_view("combination", policy.team_size(), 5);
//  Kokkos::View<int** > sgraph_distance_signature("dist sig", policy.team_size(), orbits.distance.extent(0));
//  Kokkos::View<int**> sgraph_degree_signature("Degree signature", policy.team_size(), 5);
//  Kokkos::View<int*[5][5]> induced_subgraph("Subgraph", policy.team_size());
//  Kokkos::View<bool** > visited("BFS visited", policy.team_size(), 5);
//  Kokkos::View<int**> queue("BFS queue", policy.team_size(), 5);
//  Kokkos::View<int**> distance("BFS distance", policy.team_size(), 5);
////      Kokkos::View<int[1]> node_counter("Node counter");
////      node_counter(0) = node_name;
//      Kokkos::parallel_for("Calculate GDV", policy, KOKKOS_LAMBDA(member_type team_member) {
////        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, end-node_name), [=] (int n) {
//          int node = node_name + team_member.league_rank()*team_member.team_size() + team_member.team_rank();
////          int node = Kokkos::atomic_fetch_add(&node_counter(0), 1);
////          int node = node_name + n;
////printf("Node: %d\tThreadid: %d\tTeam size: %d\tLeague size: %d\n", node, team_member.team_rank(), team_member.team_size(), team_member.league_size());
//          auto neighbor_subview = Kokkos::subview(all_neighbors, team_member.team_rank(), Kokkos::ALL());
//          auto indices_subview = Kokkos::subview(indices, team_member.team_rank(), Kokkos::ALL());
//          auto combination_subview = Kokkos::subview(combination_view, team_member.team_rank(), Kokkos::ALL());
//          auto sgraph_dist_subview = Kokkos::subview(sgraph_distance_signature, team_member.team_rank(), Kokkos::ALL());
//          auto sgraph_deg_subview = Kokkos::subview(sgraph_degree_signature, team_member.team_rank(), Kokkos::ALL());
//          auto subgraph_subview = Kokkos::subview(induced_subgraph, team_member.team_rank(), Kokkos::ALL(), Kokkos::ALL());
//          auto visited_subview = Kokkos::subview(visited, team_member.team_rank(), Kokkos::ALL());
//          auto queue_subview = Kokkos::subview(queue, team_member.team_rank(), Kokkos::ALL());
//          auto distance_subview = Kokkos::subview(distance, team_member.team_rank(), Kokkos::ALL());
//          kokkos_calculate_GDV(team_member, node, graph, orbits, neighbor_subview, indices_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, combination_counter, metrics_sa);
////        });
//      });
//      Kokkos::Experimental::contribute(metrics, metrics_sa);
//      Kokkos::View<int*> gdv_array("Send buffer for GDV", GDV_LENGTH+1);
//      for(int idx=0; idx<nrows; idx++) {
//        for(int j=0; j<GDV_LENGTH; j++) {
//          gdv_array(j) = metrics(idx, j);
//        }
//        if(idx >= nrows-1) {
//          gdv_array(GDV_LENGTH) = idx + nrows;
//        } else {
//          gdv_array(GDV_LENGTH) = idx;
//        }
////cout << "Send buffer for node " << idx << ": ";
////for(int k=0; k<gdv_array.extent(0); k++) {
////  cout << gdv_array(k) << " ";
////}
////cout << endl;
//        MPI_Send(gdv_array.data(), GDV_LENGTH+1, MPI_INT, 0, tag, MPI_COMM_WORLD);
//      }
//    } while (node_name != -1); // Exit loop if kill value is sent
//  }
//
//  //vec_calc_post_gather = MPI_Wtime() - vec_calc_start + vec_calc_post_gather;
//  vec_calc_communication_time[graph_counter - 1] = process_ends_communication - vec_calc_start;
////  cout << "Communication time on process " << rankn << " for graph " << graph_counter << " = " << vec_calc_communication_time[graph_counter - 1] << endl;
//  //vec_calc_computation_time = vec_calc_computation_end - vec_calc_computation_start + vec_calc_computation_time;
//
//  #ifdef DEBUG
//    cout << "Finished GDV Vector Calc on Rank: " << rankn << endl;
//  #endif
//
//}

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
  label += "-Rank";
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
  uint64_t k_interval = 200000000;
//  uint64_t k_interval = 4000000;

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
#ifdef INCREMENTAL
  uint32_t max_neighbors=0;
  for(int id=0; id<num_neighbors.extent(0); id++) {
    if(max_neighbors < num_neighbors(id))
      max_neighbors = num_neighbors(id);
  }
  MPI_Allreduce(MPI_IN_PLACE, &max_neighbors, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if(comm_size == 1)
  {
    Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Dynamic>> policy(1, NUM_THREADS);
    using member_type = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Dynamic>>::member_type;

#ifdef DEBUG
    cout << "Team size: " << policy.team_size() << endl;
#endif
    Kokkos::View<int**> all_neighbors("Neighbor scratch", policy.team_size(), graph.numRows());
    Kokkos::Experimental::ScatterView<uint32_t**> metrics_sa(graph_GDV);

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
#ifdef INCREMENTAL
    int capacity = max_neighbors*GDV_LENGTH;
    Kokkos::UnorderedMap<std::pair<int,int>, int> prev_map(capacity);
    Kokkos::UnorderedMap<std::pair<int,int>, int> chkpt1_map(capacity);
#endif
    for(i; i<starts.extent(0)-1; i++) {

#ifdef DIRTY_PAGE_TRACKING
      printf("Chunk %d\n", i);
      reset_dirty_bit(pid);
#endif

#ifdef AUTO_CHECKPOINT
      KokkosResilience::checkpoint(*ctx, graph_name, i, [=] () mutable {
#endif
#ifdef STANDARD
      Kokkos::View<uint32_t**> chkpt1_view("Checkpoint 1", graph_GDV.extent(0), graph_GDV.extent(1));
      Kokkos::View<uint32_t**> chkpt2_view("Checkpoint 1", graph_GDV.extent(0), graph_GDV.extent(1));
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

#ifdef STANDARD
      Kokkos::deep_copy(chkpt1_view, graph_GDV);
      std::string chkpt1_file = "checkpoint_dir/" + label + "-Checkpoint1-" + std::to_string(i);
      std::fstream chk1;
      chk1.open(chkpt1_file, std::fstream::out|std::fstream::binary);
      chk1.write((const char*)(chkpt1_view.data()), chkpt1_view.span()*sizeof(int));
      chk1.close();
      Kokkos::deep_copy(chkpt2_view, graph_GDV);
      std::string chkpt2_file = "checkpoint_dir/" + label + "-Checkpoint2-" + std::to_string(i);
      std::fstream chk2;
      chk2.open(chkpt2_file, std::fstream::out|std::fstream::binary);
      chk2.write((const char*)(chkpt2_view.data()), chkpt2_view.span()*sizeof(int));
      chk2.close();
#endif
#ifdef INCREMENTAL
      Kokkos::MDRangePolicy<Kokkos::Rank<2>> md_policy({0,0},{graph.numRows(), GDV_LENGTH});
      Kokkos::parallel_for("Checkpoint map1", md_policy, KOKKOS_LAMBDA(const int node, const int orbit) {
        std::pair<int,int> key = std::pair<int,int>(node, orbit);
        if(!prev_map.exists(key) || prev_map.value_at(prev_map.find(key)) != graph_GDV(node,orbit)) {
          auto result1 = chkpt1_map.insert(key, graph_GDV(node, orbit));
        }
      });
      std::string chkpt1_file = "incremental_checkpoint_dir/" + label + "-Checkpoint1-" + std::to_string(i);
      std::fstream chk1;
      chk1.open(chkpt1_file, std::fstream::out|std::fstream::binary);
      for(int entry=0; entry<chkpt1_map.capacity(); entry++) {
        if(chkpt1_map.valid_at(entry)) {
          auto chkpt1_key = chkpt1_map.key_at(entry);
          auto chkpt1_val = chkpt1_map.value_at(entry);
          chk1 << chkpt1_key.first << chkpt1_key.second << chkpt1_val;
        }
      }
      chk1.close();
      prev_map.clear();
      Kokkos::parallel_for(chkpt1_map.capacity(), KOKKOS_LAMBDA(uint32_t i) {
        if(chkpt1_map.valid_at(i)) {
          auto old_key = chkpt1_map.key_at(i);
          auto old_val = chkpt1_map.value_at(i);
          prev_map.insert(old_key, old_val);
        }
      });
      std::string chkpt2_file = "incremental_checkpoint_dir/" + label + "-Checkpoint2-" + std::to_string(i);
      std::fstream chk2;
      chk2.open(chkpt2_file, std::fstream::out|std::fstream::binary);
      for(int entry=0; entry<chkpt1_map.capacity(); entry++) {
        if(chkpt1_map.valid_at(entry)) {
          auto chkpt1_key = chkpt1_map.key_at(entry);
          auto chkpt1_val = chkpt1_map.value_at(entry);
          chk2 << chkpt1_key.first << chkpt1_key.second << chkpt1_val;
        }
      }
      chk2.close();
      chkpt1_map.clear();
#endif
      chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
      chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2-t1);
#ifdef DEBUG
      printf("Done with chunk: %d, time: %f\n", i, time_span.count());
#endif

#ifdef AUTO_CHECKPOINT
}, filt);
#endif
    }
  }
  else
  {
    Kokkos::deep_copy(graph_GDV, 0);
    uint32_t intervals_per_rank = num_intervals/comm_size;
    uint32_t start_index = intervals_per_rank*rankn;
    if(rankn == comm_size-1) {
      intervals_per_rank += num_intervals - (comm_size*intervals_per_rank);
    }
#ifdef DEBUG
printf("Rank %d start index: %d, intervals per rank: %d\n", rankn, start_index, intervals_per_rank);
#endif
    int offset = 0;
    Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(1, NUM_THREADS);
    using member_type = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;

    Kokkos::Experimental::ScatterView<uint32_t**> metrics_sa(graph_GDV);
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


#ifdef STANDARD 
    Kokkos::View<uint32_t**> chkpt1_view("Checkpoint 1", graph_GDV.extent(0), graph_GDV.extent(1));
    Kokkos::View<uint32_t**> chkpt2_view("Checkpoint 2", graph_GDV.extent(0), graph_GDV.extent(1));
#endif 
#ifdef INCREMENTAL
    int capacity = graph.numRows()*15;
    printf("Map capacity: %d\n", capacity);
    Kokkos::UnorderedMap<std::pair<int,int>, int> prev_map(capacity);
    Kokkos::UnorderedMap<std::pair<int,int>, int> chkpt1_map(capacity);
    Kokkos::View<uint32_t**> write_buffer("Write buffer", 3, capacity);
#endif
#ifdef HASH_DETECT
    size_t blocksize = 65536;
    size_t blocksize_elements = blocksize/4;
    Kokkos::View<uint64_t*> old_hashes("Old hashes", graph_GDV.span()/(blocksize_elements));
    Kokkos::View<uint64_t*> new_hashes("New hashes", graph_GDV.span()/(blocksize_elements));
    Kokkos::View<uint64_t*> hash_buffer("Hash buffer", graph_GDV.span()/(blocksize_elements));
    Kokkos::View<uint32_t*> write_buffer("Write buffer", graph_GDV.span());
    Kokkos::deep_copy(old_hashes, 0);
    Kokkos::deep_copy(new_hashes, 0);
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
#ifdef STANDARD
      Kokkos::deep_copy(chkpt1_view, graph_GDV);
      std::string chkpt1_file = "checkpoint_dir/" + label + "-Checkpoint1-" + std::to_string(offset);
      std::fstream chk1;
      chk1.open(chkpt1_file, std::fstream::out|std::fstream::binary);
      chk1.write((const char*)(chkpt1_view.data()), chkpt1_view.span()*sizeof(int));
      chk1.close();
      Kokkos::deep_copy(chkpt2_view, graph_GDV);
      std::string chkpt2_file = "checkpoint_dir/" + label + "-Checkpoint2-" + std::to_string(offset);
      std::fstream chk2;
      chk2.open(chkpt2_file, std::fstream::out|std::fstream::binary);
      chk2.write((const char*)(chkpt2_view.data()), chkpt2_view.span()*sizeof(int));
      chk2.close();
#endif
#ifdef INCREMENTAL
      Kokkos::View<uint64_t> buffer_counter("Buffer counter");
      buffer_counter() = 0;
      Kokkos::View<uint64_t> item_counter("Item counter");
      item_counter() = 0;
      Kokkos::View<uint64_t> repeat_counter("Repeat counter");
      repeat_counter() = 0;
printf("Starting checkpoint setup\n");
//      uint64_t b_counter = 0;
//      uint64_t non_zero = 0;
//      for(size_t _node=0; _node<graph.numRows(); _node++) {
//        for(size_t _orbit=0;_orbit<GDV_LENGTH; _orbit++) {
//          std::pair<int,int> key = std::pair<int,int>(_node, _orbit);
//          if(graph_GDV(_node,_orbit) != 0) {
//            non_zero++;
//            if(!prev_map.exists(key)) {
//              b_counter++;
//          auto result1 = chkpt1_map.insert(key, graph_GDV(_node, _orbit));
//            } else if(prev_map.value_at(prev_map.find(key)) != graph_GDV(_node,_orbit)) {
//              b_counter++;
//          auto result1 = chkpt1_map.insert(key, graph_GDV(_node, _orbit));
//            }
//          }
////          if((graph_GDV(_node,_orbit) != 0) && (!prev_map.exists(key) || (prev_map.value_at(prev_map.find(key)) != graph_GDV(_node,_orbit)))) {
////            b_counter++;
////          }
//        }
//      }
//printf("# non zeros: %lu, # of items: %lu\n", non_zero, b_counter);
//      b_counter = 0;
//      non_zero = 0;

//Kokkos::View<bool> resize("Resize");
printf("Previous map: %u\n", prev_map.size());
      Kokkos::MDRangePolicy<Kokkos::Rank<2>> md_policy({0,0},{graph.numRows(), GDV_LENGTH});
      Kokkos::parallel_for("Checkpoint map1", md_policy, KOKKOS_LAMBDA(const int node, const int orbit) {
        std::pair<int,int> key = std::pair<int,int>(node, orbit);
        if(graph_GDV(node,orbit) != 0) {
          if((!prev_map.exists(key)) || (prev_map.value_at(prev_map.find(key)) != graph_GDV(node,orbit))) {
            auto result = chkpt1_map.insert(key, graph_GDV(node, orbit));
Kokkos::atomic_add(&item_counter(), static_cast<uint64_t>(1));
//if(result.failed())
//  resize() = true;
//            uint64_t b_index = Kokkos::atomic_fetch_add(&buffer_counter(), 1);
//            write_buffer(0, b_index) = node;
//            write_buffer(1, b_index) = orbit;
//            write_buffer(2, b_index) = graph_GDV(node, orbit);
          }
          if(prev_map.exists(key) && (prev_map.value_at(prev_map.find(key)) != graph_GDV(node,orbit))) {
            Kokkos::atomic_add(&repeat_counter(), static_cast<uint64_t>(1));
          }
          uint64_t b_index = Kokkos::atomic_fetch_add(&buffer_counter(), 1);
          write_buffer(0, b_index) = node;
          write_buffer(1, b_index) = orbit;
          write_buffer(2, b_index) = graph_GDV(node, orbit);
        }
      });
printf("Rank: %d, Setup checkpoint: %lu nonzeros, %lu items, %lu repeats\n", rankn, buffer_counter(), item_counter(), repeat_counter());
      std::string chkpt1_file = "incremental_checkpoint_dir/" + label + "-Checkpoint1-" + std::to_string(offset);
      std::fstream chk1;
      chk1.open(chkpt1_file, std::fstream::out|std::fstream::binary);
      chk1.write((const char*)(write_buffer.data()), buffer_counter()*3*sizeof(int));
//      for(int entry=0; entry<chkpt1_map.capacity(); entry++) {
//        if(chkpt1_map.valid_at(entry)) {
//          auto chkpt1_key = chkpt1_map.key_at(entry);
//          auto chkpt1_val = chkpt1_map.value_at(entry);
//          chk1 << chkpt1_key.first << chkpt1_key.second << chkpt1_val;
//        }
//      }
      chk1.close();
//printf("Capacity: %d, size: %d\n", chkpt1_map.capacity(), chkpt1_map.size());
//      std::string chkpt2_file = "incremental_checkpoint_dir/" + label + "-Checkpoint2-" + std::to_string(offset);
//      std::fstream chk2;
//      chk2.open(chkpt2_file, std::fstream::out|std::fstream::binary);
//      chk2.write((const char*)(write_buffer.data()), buffer_counter()*3*sizeof(int));
//      for(int entry=0; entry<chkpt1_map.capacity(); entry++) {
//        if(chkpt1_map.valid_at(entry)) {
//          auto chkpt1_key = chkpt1_map.key_at(entry);
//          auto chkpt1_val = chkpt1_map.value_at(entry);
//          chk2 << chkpt1_key.first << chkpt1_key.second << chkpt1_val;
//        }
//      }
//      chk2.close();
      Kokkos::deep_copy(prev_map, chkpt1_map);
//      prev_map.clear();
//      Kokkos::parallel_for(chkpt1_map.capacity(), KOKKOS_LAMBDA(uint32_t i) {
//        if(chkpt1_map.valid_at(i)) {
//          auto old_key = chkpt1_map.key_at(i);
//          auto old_val = chkpt1_map.value_at(i);
//          prev_map.insert(old_key, old_val);
//        }
//      });
      chkpt1_map.clear();
#endif
#ifdef HASH_DETECT
      generate_hashes(graph_GDV, new_hashes, blocksize);
//      chk1 << new_hashes.span();
//      chk1.write((const char*)(new_hashes.data()), new_hashes.span()*sizeof(uint64_t));
      Kokkos::View<uint64_t> buffer_counter("Buffer counter");
      buffer_counter() = 0;
      Kokkos::parallel_for("Find dirty pages", Kokkos::RangePolicy<>(0, new_hashes.span()), KOKKOS_LAMBDA(const size_t idx) {
//      for(size_t idx=0; idx<new_hashes.span(); idx++) {
        if(old_hashes(idx) != new_hashes(idx)) {
//          chk1 << new_hashes(idx);
//          chk1.write((const char*)(graph_GDV.data()+(idx*blocksize_elements)), blocksize);
          uint32_t b_index = Kokkos::atomic_fetch_add(&buffer_counter(), 1);
          hash_buffer(b_index) = new_hashes(idx);
          std::memcpy((write_buffer.data()+(b_index*blocksize_elements)), (graph_GDV.data()+idx*blocksize_elements), blocksize);
//          for(int i=0; i<blocksize_elements; i++) {
//            write_buffer(b_index*blocksize_elements + i) = *(graph_GDV.data() + idx*blocksize_elements + i);
//          }
        } else {
          printf("Clean page %d\n", idx);
        }
//      }
      });
//if(offset > 0) {
      std::string chkpt1_file = "incremental_checkpoint_dir/" + label + "-Checkpoint1-" + std::to_string(offset);
      std::fstream chk1;
      chk1.open(chkpt1_file, std::fstream::out|std::fstream::binary);
      chk1 << new_hashes.span();
      chk1.write((const char*)(new_hashes.data()), new_hashes.span()*sizeof(uint64_t));
      chk1 << buffer_counter();
      chk1.write((const char*)(hash_buffer.data()), buffer_counter()*sizeof(uint64_t));
      chk1.write((const char*)(write_buffer.data()), buffer_counter()*blocksize_elements*sizeof(uint32_t));
      chk1.close();
//}
      Kokkos::deep_copy(old_hashes, new_hashes);
      Kokkos::deep_copy(new_hashes, 0);
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



//  else if (rankn == 0) 
//  {
//    int i;
//    if (num_intervals < comm_size) {
//      // Send all nodes if comm size is bigger
//      int range[2];
//      for (i = 0; i < starts.extent(0)-1; i++) {
//        range[0] = i;
//        range[1] = i+CHUNK_SIZE;
//        MPI_Send(range, 2, MPI_INT, i+1, tag, MPI_COMM_WORLD);
//      }
//
//      // Send termination to finish processes
//      int flag;
//      for (i = 1; i < comm_size; i++) {
////        flag = -1;
////        MPI_Send(&flag, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
//        range[0] = -1;
//        range[1] = -1;
//        MPI_Send(range, 2, MPI_INT, i, tag, MPI_COMM_WORLD);
//      }
//
//    } else { // There are more nodes in graph than there are MPI processes
//
//      // First get each process busy
//      int send_node; // Corresponds to index of node to send
//      int rcv_node;  // Corresponds to name of graph node recieved from worker rank
//      int range[2];
//      int rec_count = 0;
//
//#ifdef AUTO_CHECKPOINT
//      auto ctx = KokkosResilience::make_context(MPI_COMM_SELF, "/home/ntan1/Src_Fido_Kokkos/fido.json");
//      printf("Created root context\n");
//      const auto filt = KokkosResilience::Filter::NthIterationFilter(1);
//      printf("Created root filter\n");
//      rec_count = KokkosResilience::latest_version(*ctx, label);
//      if(rec_count < 0)
//        rec_count = 0;
//      printf("Got latest counter %d\n", rec_count);
//#endif
//
//      for (i=1; i < comm_size; i++) {
//	      send_node = (i-1)*CHUNK_SIZE;
//        #ifdef DEBUG
//          cout << "Sending node " << send_node << " to rank " << i << endl;
//        #endif
////        MPI_Send(&send_node, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
//        range[0] = send_node;
//        range[1] = send_node+CHUNK_SIZE;
//        MPI_Send(range, 2, MPI_INT, i, tag, MPI_COMM_WORLD);
//      }
//
////cout << "Sent initial batch of nodes\n";
//
//      // Start probing and recieving results from other processes
//      send_node = (comm_size-1)*CHUNK_SIZE;
//      do {
//	
//        // First probe for completed work
//      	int flag;
//      	MPI_Status master_status;
//      	MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &master_status);
//      
//      	if (flag == 1) {
//      	
//      	  // Recieve gdv from finished process
//      	  i = master_status.MPI_SOURCE;
//          Kokkos::View<int[GDV_LENGTH+1]> gdv_array_d("Recv buffer");
//          auto gdv_array_h = Kokkos::create_mirror_view(gdv_array_d);
//      	  MPI_Recv(gdv_array_h.data(), GDV_LENGTH + 1, MPI_INT, i, tag, MPI_COMM_WORLD, &master_status);
//          Kokkos::deep_copy(gdv_array_d, gdv_array_h);
//
//      	  #ifdef DEBUG
//            cout << "Recieved GDV for node " << rcv_node << " from rank " << i << ": " << endl;
//            cout << gdv_array_h(GDV_LENGTH) << ": ";
//            for (int j = 0; j < GDV_LENGTH; j++) {
//              cout << gdv_array_h(j) << ", ";
//            }
//            cout << endl;
//          #endif
//      
//      	  // Organize recieved data into returning array
//          rcv_node = gdv_array_h(GDV_LENGTH);
//          // We are updating each vertex in the nodes neighborhood.
//          // Last GDV update for the node will be node+graph_size 
//          // as a signal for when we're done with this node
//          if(rcv_node < 0) {
//            rcv_node *= -1;
//            rec_count += CHUNK_SIZE;
//            int range[2];
////#ifdef AUTO_CHECKPOINT
////KokkosResilience::checkpoint(*ctx, label, rec_count, [=] () mutable {
////#endif
//            if(send_node < num_intervals) {
//              range[0] = send_node;
//              range[1] = send_node+CHUNK_SIZE;
//      	      MPI_Send(range, 2, MPI_INT, i, tag, MPI_COMM_WORLD);
//      	      send_node += CHUNK_SIZE;
//      	    } else { // Send termination
//              range[0] = -1;
//              range[1] = -1;
//printf("Sent kill signal to %d\n", i);
//              MPI_Send(range, 2, MPI_INT, i, tag, MPI_COMM_WORLD);
//            }
//            Kokkos::parallel_for("Update GDV", Kokkos::RangePolicy<>(0, GDV_LENGTH), KOKKOS_LAMBDA(const int j) {
//              graph_GDV(rcv_node, j) += gdv_array_d(j);
//            });
////#ifdef AUTO_CHECKPOINT
////}, filt);
////#endif
//          } else {
//            Kokkos::parallel_for("Update GDV", Kokkos::RangePolicy<>(0, GDV_LENGTH), KOKKOS_LAMBDA(const int j) {
//              graph_GDV(rcv_node, j) += gdv_array_d(j);
//            });
//          }
//      	}
//      } while (rec_count < num_intervals);
//printf("Received all GDVs\n");
//
//      process_ends_communication = MPI_Wtime();
//      //vec_calc_prior_gather = MPI_Wtime() - vec_calc_start + vec_calc_prior_gather;
//
//      // Sort return vector
////      sort(graph_GDV->begin(), graph_GDV->end(), node_id_order());
//      #ifdef DEBUG
//        cout << "Constructed return GDV array" << endl;
//        for (i = 0; i < graph_GDV->size(); i++) {
//      	  cout << graph_GDV->at(i).node << ": ";
//      	  for (int j = 0; j < graph_GDV->at(i).GDV.size(); j++) {
//      	    cout << graph_GDV->at(i).GDV[j] << ", ";
//      	  }
//      	  cout << endl;
//      	}	
//      #endif
//    }
//printf("Rank 0 done\n");
//  }
//  else // Instructions for work processes
//  { 
//    int node_name;
//    MPI_Status worker_status;
//    int chunk_idx;
//    int chunk_counter = 0;
//#ifdef AUTO_CHECKPOINT
//    auto ctx = KokkosResilience::make_context(MPI_COMM_SELF, "/home/ntan1/Src_Fido_Kokkos/fido.json");
//    printf("Created context\n");
//    const auto filt = KokkosResilience::Filter::NthIterationFilter(1);
//    printf("Created filter\n");
//    chunk_counter = KokkosResilience::latest_version(*ctx, label);
//    if(chunk_counter < 0)
//      chunk_counter = 0;
//    printf("Got latest counter %d\n", chunk_idx);
//#endif
//
//    Kokkos::View<int[2]> range("Start,end");
//    do {
//
////cout << "Starting worker thread\n";
////      int range[2];
//      MPI_Recv(range.data(), 2, MPI_INT, 0, tag, MPI_COMM_WORLD, &worker_status);
//      node_name = range(0);
////      MPI_Recv(&node_name, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &worker_status);
////cout << "Got new batch of nodes starting at " << node_name << "\n";
//      #ifdef DEBUG
//        cout << "Recieved node " << node_name << " at rank " << rankn << endl;
//      #endif
//      if (node_name == -1) {
////  MPI_Comm_split(MPI_COMM_WORLD, 0, rankn, &worker_comm);
////printf("Split communicator\n");
//        process_ends_communication = MPI_Wtime();
//        //vec_calc_prior_gather = MPI_Wtime() - vec_calc_start + vec_calc_prior_gather;
////cout << "Thread " << rankn << " done\n";
//        break;
//      } 
////      else 
////      {
////  MPI_Comm_split(MPI_COMM_WORLD, 1, rankn, &worker_comm);
////printf("Split communicator\n");
////      }
////      int end = node_name+CHUNK_SIZE;
//      int end = range(1);
//      int nrows = graph.numRows();
////      if(end > nrows)
////        end = nrows;
//      if(end > num_intervals)
//        end = num_intervals;
//printf("Received chunk (%d,%d)\n", node_name, end);
//#ifdef INSTRUMENT
////printf("Combination | # of updated regions | avg size of region | (# of identical regions, total size)\n");
////printf("Iteration|# of regions|Total size|Avg size of region|%% of entries updated|(updates,empty)|# nonzero|# contiguous regions\n");
//#endif
//
//      Kokkos::deep_copy(graph_GDV, 0);
//      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(1, NUM_THREADS);
//      using member_type = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
//
////cout << "Team size: " << policy.team_size() << endl;
//      Kokkos::Experimental::ScatterView<int**> metrics_sa(graph_GDV);
//
//      Kokkos::View<int*> combination_counter("Per thread counter", policy.team_size());
//      Kokkos::deep_copy(combination_counter, 0);
//      Kokkos::View<int**> all_neighbors("Neighbor scratch", policy.team_size(), graph.numRows());
//      Kokkos::View<int**> indices("Index", policy.team_size(), 5);
//      Kokkos::View<int**> combination_view("combination", policy.team_size(), 5);
//      Kokkos::View<int**> sgraph_distance_signature("dist sig", policy.team_size(), orbits.distance.extent(0));
//      Kokkos::View<int**> sgraph_degree_signature("Degree signature", policy.team_size(), 5);
//      Kokkos::View<int*[5][5]> induced_subgraph("Subgraph", policy.team_size());
//      Kokkos::View<bool** > visited("BFS visited", policy.team_size(), 5);
//      Kokkos::View<int**> queue("BFS queue", policy.team_size(), 5);
//      Kokkos::View<int**> distance("BFS distance", policy.team_size(), 5);
////    int chunk_idx;
////#ifdef AUTO_CHECKPOINT
//////    auto ctx = KokkosResilience::make_context(worker_comm, "/home/ntan1/Src_Fido_Kokkos/fido.json");
////    auto ctx = KokkosResilience::make_context(MPI_COMM_SELF, "/home/ntan1/Src_Fido_Kokkos/fido.json");
////    printf("Created context\n");
////    const auto filt = KokkosResilience::Filter::NthIterationFilter(1);
////    printf("Created filter\n");
////    chunk_idx = KokkosResilience::latest_version(*ctx, graph_name);
////    if(chunk_idx < 0)
////      chunk_idx = 0;
////    printf("Got latest counter %d\n", chunk_idx);
////#endif
//      for(int chunk_idx=node_name; chunk_idx<end; chunk_idx++) {
//chunk_counter += 1;
//#ifdef AUTO_CHECKPOINT
//KokkosResilience::checkpoint(*ctx, label, chunk_counter, [=] () mutable {
//#endif
//printf("Rank %d: chunk counter: %d\n", rankn, chunk_counter-1);
////printf("Node group %d: (%d,%d)\n", chunk_idx, starts(chunk_idx), starts(chunk_idx+1));
//        bool multi_node = (chunk_idx+1 == starts.extent(0) || starts(chunk_idx+1) != starts(chunk_idx)) && 
//                          (chunk_idx==0 || starts(chunk_idx-1) != starts(chunk_idx));
//        if(multi_node) {
//          int start_node = starts(chunk_idx);
//          int end_node;
//          if(chunk_idx+1 == starts.extent(0)) {
//            end_node = graph.numRows();
//          } else {
//            end_node = starts(chunk_idx+1);
//          }
//          Kokkos::parallel_for("Calculate GDV", policy, KOKKOS_LAMBDA(member_type team_member) {
//            Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, end_node-start_node), [=] (int n) {
//              int node = start_node + n;
//              auto neighbor_subview = Kokkos::subview(all_neighbors, team_member.team_rank(), Kokkos::ALL());
//              auto indices_subview = Kokkos::subview(indices, team_member.team_rank(), Kokkos::ALL());
//              auto combination_subview = Kokkos::subview(combination_view, team_member.team_rank(), Kokkos::ALL());
//              auto sgraph_dist_subview = Kokkos::subview(sgraph_distance_signature, 
//                                                          team_member.team_rank(), 
//                                                          Kokkos::ALL());
//              auto sgraph_deg_subview = Kokkos::subview(sgraph_degree_signature, 
//                                                        team_member.team_rank(), 
//                                                        Kokkos::ALL());
//              auto subgraph_subview = Kokkos::subview(induced_subgraph, 
//                                                      team_member.team_rank(), 
//                                                      Kokkos::ALL(), 
//                                                      Kokkos::ALL());
//              auto visited_subview = Kokkos::subview(visited, team_member.team_rank(), Kokkos::ALL());
//              auto queue_subview = Kokkos::subview(queue, team_member.team_rank(), Kokkos::ALL());
//              auto distance_subview = Kokkos::subview(distance, team_member.team_rank(), Kokkos::ALL());
//              kokkos_calculate_GDV(team_member, node, graph, orbits, neighbor_subview, indices_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, combination_counter, metrics_sa);
//            });
//          });
//          Kokkos::Experimental::contribute(graph_GDV, metrics_sa);
//          metrics_sa.reset();
//        } else {
//          auto neighbor_scratch = Kokkos::subview(all_neighbors, 0, Kokkos::ALL());
//          int n_neighbors = EssensKokkos::find_neighbours(starts(chunk_idx), graph, 4, neighbor_scratch);
//          auto neighbor_subview = Kokkos::subview(all_neighbors, 0, std::pair<int,int>(0,n_neighbors));
////printf("Number of neighbors: %d\n", n_neighbors);
////for(int n = 0; n<n_neighbors; n++) {
////  printf("%d ", neighbor_subview(n));
////}
////printf("\n");
//          int start_combination;
//          int end_combination;
//          int start_comb_subgraph[5] = {0,0,0,0,0};
//          int end_comb_subgraph[5] = {0,0,0,0,0};
//          start_comb_subgraph[0] = 0;
//          end_comb_subgraph[0] = 0;
//          int start_chunk = 0;
//          for(int j=0; j<chunk_idx; j++) {
//            if(starts(j) == starts(chunk_idx))
//              start_chunk++;
//          }
//          bool first_chunk = (chunk_idx == 0) || 
//                              starts(chunk_idx-1) != starts(chunk_idx) && 
//                              starts(chunk_idx+1) == starts(chunk_idx);
//          bool last_chunk = starts(chunk_idx-1) == starts(chunk_idx) && 
//                            starts(chunk_idx+1) != starts(chunk_idx);
//          bool middle_chunk = starts(chunk_idx-1) == starts(chunk_idx) && 
//                              starts(chunk_idx+1) == starts(chunk_idx);
////printf("Start chunk: %d\n", start_chunk);
//          if(last_chunk) {
//            // Last chunk
////printf("Last chunk of %d\n", starts(chunk_idx));
//            start_combination = start_chunk*k_interval;
//            end_combination = num_combinations(starts(chunk_idx));
//            int64_t counter = end_combination;
//            for(int j=4; j>0; j--) {
//              int64_t n_comb = get_num_combinations(n_neighbors, j);
//              end_comb_subgraph[j] = n_comb;
//              counter -= n_comb;
//              if(counter > start_combination) {
//                start_comb_subgraph[j] = 0;
//              } else {
//                start_comb_subgraph[j] = start_combination-counter;
//                break;
//              }
//            }
//          } else if(first_chunk) {
//            // First chunk
////printf("First chunk of %d\n", starts(chunk_idx));
//            start_combination = 0;
//            end_combination = k_interval;
//            size_t counter = k_interval;
//            for(int j=1; j<5; j++) {
//              int64_t n_comb = get_num_combinations(n_neighbors, j);
//              if(counter > n_comb) {
//                end_comb_subgraph[j] = n_comb;
//                counter -= n_comb;
//              } else {
//                end_comb_subgraph[j] = counter;
//                break;
//              }
//            }
//          } else if(middle_chunk) {
//            // Middle chunk
////printf("Middle chunk of %d\n", starts(chunk_idx));
//            start_combination = start_chunk*k_interval;
//            end_combination = start_combination+k_interval;
//            size_t counter = 0;
//            for(int j=1; j<5; j++) {
//              int64_t n_comb = get_num_combinations(n_neighbors, j);
//              if(start_combination > counter && start_combination < counter+n_comb) {
//                start_comb_subgraph[j] = start_combination-counter;
//              }
//              if(end_combination > counter+n_comb && counter+n_comb >= start_combination) {
//                end_comb_subgraph[j] = n_comb;
//              }
//              if(end_combination > counter && end_combination < counter+n_comb) {
//                end_comb_subgraph[j] = end_combination-counter;
//                break;
//              }
//              counter += n_comb;
//            }
//          }
//          
////printf("Start combination: %d\n", start_combination);
////printf("End combination: %d\n", end_combination);
////for(int idx=0; idx<5; idx++) {
////  printf("%d ", start_comb_subgraph[idx]);
////}
////printf("\n");
////for(int idx=0; idx<5; idx++) {
////  printf("%d ", end_comb_subgraph[idx]);
////}
////printf("\n");
//          
//          int node = starts(chunk_idx);
//          Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> team_policy(1, NUM_THREADS);
//          Kokkos::View<int**> scratch_block("Indices scratch", team_policy.team_size(), graph.numRows());
//          for(int node_count = 1; node_count < 5; node_count++) {
////printf("Node count: %d\n", node_count);
//            int64_t s_comb = start_comb_subgraph[node_count];
//            int64_t e_comb = end_comb_subgraph[node_count];
//            if(e_comb-s_comb > 0) {
////printf("Node %d: subgraph size: %d, start,end combination (%d, %d)\n", node, node_count, s_comb, e_comb);
//              Kokkos::parallel_for("Calculate GDV", team_policy, KOKKOS_LAMBDA(member_type team_member) {
//                Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, e_comb-s_comb), [=] (int64_t idx) {
//                  int64_t combination_num = idx+s_comb;
////if(idx == 0 || idx == e_comb-s_comb-1)
////  printf("Running combination %ld\n", combination_num);
//                  auto indices_subview = Kokkos::subview(indices, team_member.team_rank(), Kokkos::ALL());
//                  auto combination_subview = Kokkos::subview(combination_view, team_member.team_rank(), Kokkos::ALL());
//                  auto sgraph_dist_subview = Kokkos::subview(sgraph_distance_signature, 
//                                                              team_member.team_rank(), 
//                                                              Kokkos::ALL());
//                  auto sgraph_deg_subview = Kokkos::subview(sgraph_degree_signature, 
//                                                            team_member.team_rank(), 
//                                                            Kokkos::ALL());
//                  auto subgraph_subview = Kokkos::subview(induced_subgraph, 
//                                                          team_member.team_rank(), 
//                                                          Kokkos::ALL(), 
//                                                          Kokkos::ALL());
//                  auto visited_subview = Kokkos::subview(visited, team_member.team_rank(), Kokkos::ALL());
//                  auto queue_subview = Kokkos::subview(queue, team_member.team_rank(), Kokkos::ALL());
//                  auto distance_subview = Kokkos::subview(distance, team_member.team_rank(), Kokkos::ALL());
//                  auto scratch_view = Kokkos::subview(scratch_block, team_member.team_rank(), Kokkos::ALL());
//                  combination_from_position(scratch_view, combination_num, n_neighbors, node_count);
//                  for(int j=0; j<node_count; j++) {
//                      combination_subview(j) = neighbor_subview(scratch_view(j));
//                  }
//                  combination_subview(node_count) = node;
//                  kokkos_calculate_GDV(team_member, node, node_count, graph, orbits, neighbor_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, metrics_sa);
//                });
//              });
//            }
//          }
//          Kokkos::Experimental::contribute(graph_GDV, metrics_sa);
//          metrics_sa.reset();
//        }
//#ifdef AUTO_CHECKPOINT
//}, filt);
//#endif
//      }
//      Kokkos::View<int*> gdv_array("Send buffer for GDV", GDV_LENGTH+1);
//      for(int idx=0; idx<nrows; idx++) {
//        for(int j=0; j<GDV_LENGTH; j++) {
//          gdv_array(j) = graph_GDV(idx, j);
//        }
//        if(idx >= nrows-1) {
//          gdv_array(GDV_LENGTH) = -idx;
//        } else {
//          gdv_array(GDV_LENGTH) = idx;
//        }
//        MPI_Send(gdv_array.data(), GDV_LENGTH+1, MPI_INT, 0, tag, MPI_COMM_WORLD);
//      }
//    } while (node_name != -1); // Exit loop if kill value is sent
//  }

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
#ifdef INSTRUMENT
  int n_updates = 0;
  int num_blank_updates = 0;
  int num_meaningful_updates = 0;
#endif
  auto gdvMetrics = gdvMetrics_sa.access();
//  int combination_count = 0;
//  int k_interval = 1000000;
  int num_neighbors = EssensKokkos::find_neighbours(node, graph, 4, neighbor_buff);
  auto neighbors = Kokkos::subview(neighbor_buff, std::pair<int,int>(0, num_neighbors));
//  printf("Allocated and found neighbors of %d: %d total neighbors\n", node, neighbors.size());
//  int iter_counter = counter(team_member.team_rank());
#ifdef INSTRUMENT
Kokkos::View<bool**> gdv_updated("Updated entries", graph.numRows(), orbits.distance.extent(0));
Kokkos::View<bool**> gdv_updated_prev("Previously updated entries", graph.numRows(), orbits.distance.extent(0));
Kokkos::deep_copy(gdv_updated, false);
Kokkos::deep_copy(gdv_updated_prev, false);
int valid_combinations = 0;
chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
#endif
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
  
#ifdef INSTRUMENT
n_updates+=1;
gdv_updated(combination(v), i) = true;
#endif
            }
          }
        }
#ifdef INSTRUMENT
valid_combinations += 1;
#endif
      }
      generator.kokkos_next(indices);
//if(node == 1) {
//if(comb_counter < 100) {
//  printf("Number of updates for iter %d: %d\n", comb_counter, n_comb);
//  n_comb = 0;
//  comb_counter++;
//}
//}

#ifdef INSTRUMENT
if(n_updates > 0) {
  bool contiguous = false;
  int start_node = 0;
  int start_orbit = 0;
  int num_big_regions = 0;
  int total_region_size = 0;
  int num_regions = 0;
  int num_updates = 0;
  int total_updates = 0;
  int new_updates = 0;
  for(int i=0; i<gdv_updated.extent(0); i++) {
    for(int j=0; j<gdv_updated.extent(1); j++) {
      if(gdv_updated(i,j) || gdv_updated_prev(i,j))
        total_updates++;
      if(gdv_updated(i,j) && !gdv_updated_prev(i,j))
        new_updates += 1;
      if(gdv_updated(i,j)) {
        num_updates += 1;
      }
      if(gdv_updated(i,j) && !contiguous) {
        contiguous = true;
        total_region_size += 1;
//        num_updates += 1;
        start_node = i;
        start_orbit = j;
      } else if(gdv_updated(i,j) && contiguous) {
        total_region_size += 1;
//        num_updates += 1;
      } else if(!gdv_updated(i,j) && contiguous) {
        contiguous = false;
        num_regions += 1;
        if(j != start_orbit+1) {
          num_big_regions += 1;
        }
      }
    }
  }
  if(num_updates > 0) {
    num_meaningful_updates+=1;
  } else {
    num_blank_updates+=1;
  }
//if(combination_count % k_interval == 0 && num_updates > 0) {
if(iter_counter % k_interval == 0 && num_updates > 0) {
// Combination iteration
    printf("%5d\t", iter_counter);
// Combination
//    cout << "(";
//    for(int i=0; i<combination.size()-1; i++) {
//      printf("%4d, ", combination(i));
//    }
//    printf("%4d)", combination(combination.size()-1));
//    for(int i=0; i<5-combination.size(); i++) {
//      printf("      ");
//    }
// Number of regions. Includes regions of size 1
//    cout << "\t(" << num_regions << ")\t";
// Total # of updated values
//    printf("(%d)\t", total_region_size);
    printf("(%d)\t", total_updates);
// Average size of region
//    if(num_regions == 0) {
//      printf("(%7f)\t", 0.0);
//    } else {
//      printf("(%7f)\t", static_cast<float>(total_region_size)/static_cast<float>(num_regions));
//    }
// Percentage of GDV entries updated
    printf("%f\t", static_cast<float>(total_updates)/static_cast<float>(gdv_updated.size()));
// Number of updates
    printf("(%d)\t", num_updates);
// Percentage of updates
    printf("%f\t", static_cast<float>(num_updates)/static_cast<float>(gdv_updated.size()));
// Number of updates
    printf("(%d)\t", new_updates);
// Percentage of updates
    printf("%f\t", static_cast<float>(new_updates)/static_cast<float>(gdv_updated.size()));
//    total_region_size = 0;
//    num_regions = 0;
//    contiguous = false;
//    for(int i=0; i<gdv_updated.extent(0); i++) {
//      for(int j=0; j<gdv_updated.extent(1); j++) {
//        if(gdv_updated(i,j) && gdv_updated_prev(i,j) && !contiguous) {
//          contiguous = true;
//          total_region_size += 1;
//        } else if(gdv_updated(i,j) && gdv_updated_prev(i,j) && contiguous) {
//          total_region_size += 1;
//        } else if((!gdv_updated(i,j) || !gdv_updated_prev(i,j)) && contiguous) {
//          contiguous = false;
//          num_regions += 1;
//        }
//      }
//    }
//// # of regions and size of regions that are shared with the previous interval
//    cout << "(" << num_regions << "," << total_region_size << ")";
//    Kokkos::deep_copy(gdv_updated_prev, gdv_updated);
//    Kokkos::deep_copy(gdv_updated, 0);
//// Percentage of previously updated GDV entries
//    printf("\t%f", static_cast<float>(num_updates)/static_cast<float>(gdv_updated.size()));
// # of meaningful updated vs empty updates
//    printf("\t(%d,%d)", num_meaningful_updates, num_blank_updates);
// # of valid combinations
//    printf("\t%d", valid_combinations);
// # of regions larger than 1
//    printf("\t%d", num_big_regions);
    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2-t1);
    printf("%fs\t", time_span.count());
    printf("\n");
    for(int i=0; i<gdv_updated.extent(0); i++) {
      for(int j=0; j<gdv_updated.extent(1); j++) {
        gdv_updated_prev(i,j) = gdv_updated_prev(i,j) || gdv_updated(i,j);
        gdv_updated(i,j) = false;
      }
    }
}
}
n_updates = 0;
if(iter_counter % k_interval == 0) {
t1 = chrono::steady_clock::now();
}
#endif
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

