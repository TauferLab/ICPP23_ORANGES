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
#include <cstdlib>
#include <ctime>
#include <mpi.h>
#include <fstream>
#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include <Kokkos_ScatterView.hpp>

#define GDV_LENGTH 73
#define CHUNK_SIZE 64
#define SIM_MAT_CUTOFF SIM_MAT_CUT
//#define DEBUG
//if (TRACK_RUNTIME) {
#define RUNTIME TRACK_RUNTIME
#define RUNTIME_VAL 1
//}
#define NUM_THREADS THREAD_COUNT

/*void Calculate_GDV(int ,A_Network ,vector<OrbitMetric>&, GDVMetric&);
void Calculate_GDV(int node,A_Network Graph,vector<OrbitMetric> &orbits, vector<GDVMetric> &gdvMetrics);
void readin_orbits(  ifstream* ,vector<OrbitMetric>* );
void convert_string_vector_int(string* , vector<int>* ,string );
void Similarity_Metric_calculation_for_two_graphs(A_Network, A_Network,vector<OrbitMetric>, string, string);
double GDV_distance_calculation(GDVMetric&, GDVMetric&);
void metric_formula(GDVMetric&, double*);
void GDV_vector_calculation(A_Network,vector<GDVMetric>*,  vector<OrbitMetric>, const char*, int);*/

void kokkos_readin_orbits(ifstream *file, Orbits& orbits );
void readin_graph(ifstream* file, matrix_type& graph);
void kokkos_Similarity_Metric_calculation_for_two_graphs(const matrix_type&, const matrix_type&, const Orbits&, string, string, int);
void convert_string_vector_int(string* , vector<int>* ,string );
template <class SubviewType>
KOKKOS_INLINE_FUNCTION double kokkos_GDV_distance_calculation(SubviewType&, SubviewType&);
template <class SubviewType>
KOKKOS_INLINE_FUNCTION double kokkos_metric_formula(SubviewType &gdvm);
void kokkos_GDV_vector_calculation(const matrix_type&, GDVs&, const Orbits&, const char*, int, int);
template<class NeighborView, class IntView, class GraphView, class BoolView, class CounterView
>
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
                      CounterView& counters,
                      Kokkos::Experimental::ScatterView<int**> gdvMetrics_sa
                    );

// Functions to add
//void record_calculated_results();
void dynamic_serial_work(const matrix_type&, GDVs&, const Orbits&, Kokkos::View<unsigned long int*>, Kokkos::View<unsigned long int*>, int, int);
void dynamic_master_work(GDVs&, int, int);
void dynamic_worker_work(const matrix_type&, const Orbits&, int, int);

struct node_id_order {
  inline bool operator() (const GDVMetric& gdv_met_1, const GDVMetric& gdv_met_2) {
    return (gdv_met_1.node < gdv_met_2.node);
  }
};

// Define variables for keeping track of time for load imbalancing tests.
#if RUNTIME == RUNTIME_VAL
double total_time_taken;
double vec_calc_communication_time[2] = {};
double pre_process_time;
double report_output_time;
double* vec_calc_computation_time_X;
double* vec_calc_computation_time_Y;
int* vec_calc_thread_assign_X;
int* vec_calc_thread_assign_Y;
int* vec_calc_proc_assign_X;
int* vec_calc_proc_assign_Y;
double* chunk_gdvs_computation_time_X;
double* chunk_gdvs_computation_time_Y;
int* chunk_gdvs_proc_assign_X;
int* chunk_gdvs_proc_assign_Y;
void record_runtimes(const matrix_type&, const matrix_type&, double*, int, time_t, int, int);
#endif
//double vec_calc_prior_gather;
//double vec_calc_post_gather;
//int vec_calc_avg_node_deg;

using namespace std;

int main(int argc, char *argv[]) {

//  MPI_Init(&argc,&argv);
  int provided;
  int num_threads = atoi(argv[4]);
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided); 
  Kokkos::initialize(argc, argv);
  {
  //  if (rank == 0) {
  //clock_t out_tStart = clock();
  time_t now = time(0);
  #if RUNTIME == RUNTIME_VAL
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

//  vector<OrbitMetric> orbits;
//  readin_orbits(&the_file2,&orbits);
  Orbits k_orbits;
  kokkos_readin_orbits(&the_file2, k_orbits);

  //A_Network X;
  //A_Network Y;
  //readin_network(&X,argv[1],0,-1);
  //readin_network(&Y,argv[2],0,-1);
  GDV_functions test_gdvf;

  matrix_type graphx, graphy;
  readin_graph(&the_file0, graphx);
  readin_graph(&the_file1, graphy);

//for(int i=0; i<X.size(); i++) {
//  auto row = graphx.rowConst(i);
////cout << "A_Network: " << X[i].ListW.size() << " Graph: " << row.length << " row: " << X[i].Row << " " << i << endl;
//  for(int j=0; j<X[i].ListW.size(); j++) {
//    int xv = X[i].ListW[j].first;
//    int gxv = row.colidx(j);
//    if(xv != gxv) {
//      cout << "Graphs don't match\n";
//      cout << xv << " vs " << gxv << endl;
//      break;
//    }
//  }
//}
//for(int i=0; i<Y.size(); i++) {
//  auto row = graphy.rowConst(i);
////cout << "A_Network: " << X[i].ListW.size() << " Graph: " << row.length << " row: " << X[i].Row << " " << i << endl;
//  for(int j=0; j<Y[i].ListW.size(); j++) {
//    int xv = Y[i].ListW[j].first;
//    int gxv = row.colidx(j);
//    if(xv != gxv) {
//      cout << "Graphs don't match\n";
//      cout << xv << " vs " << gxv << endl;
//      break;
//    }
//  }
//}

  //clock_t out_tStart = clock();

  int numtasks, rank, dest, source, rc, count, tag=0;
  MPI_Status Stat;   // required variable for receive routines                                                                                                                                                                          

  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
if(rank == 0) 
  printf("Comm size: %d\n", numtasks);

  // Allocate space for time recording by graph node
#if RUNTIME == RUNTIME_VAL
  vec_calc_computation_time_X = new double[graphx.numRows()]();
  vec_calc_computation_time_Y = new double[graphy.numRows()]();
  vec_calc_thread_assign_X = new int[graphx.numRows()]();
  vec_calc_thread_assign_Y = new int[graphy.numRows()]();
  vec_calc_proc_assign_X = new int[graphx.numRows()]();
  vec_calc_proc_assign_Y = new int[graphy.numRows()]();
  chunk_gdvs_computation_time_X = new double[graphx.numRows()/CHUNK_SIZE + 1]();
  chunk_gdvs_computation_time_Y = new double[graphy.numRows()/CHUNK_SIZE + 1]();
  chunk_gdvs_proc_assign_X = new int[graphx.numRows()/CHUNK_SIZE + 1]();
  chunk_gdvs_proc_assign_Y = new int[graphy.numRows()/CHUNK_SIZE + 1]();
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
if(rank == 0) {
  cout << "# of vertices: " << graphx.numRows() << endl;
  cout << "# of vertices: " << graphy.numRows() << endl;
}

  #if RUNTIME == RUNTIME_VAL
  pre_process_time = (double)(MPI_Wtime() - pre_tStart);
  #endif
  // Perform Similarity Calculations
//  Similarity_Metric_calculation_for_two_graphs(X,Y,orbits, graph_name1, graph_name2);
  kokkos_Similarity_Metric_calculation_for_two_graphs(graphx, graphy, k_orbits, graph_name1, graph_name2, num_threads);

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
#if RUNTIME == RUNTIME_VAL
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
  #if RUNTIME == RUNTIME_VAL
    record_runtimes(graphx, graphy, time_buff, num_times, now, rank, numtasks);
  #endif

  //printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  //cout << "Time taken on rank " << rank << " = " << (double)(clock() - tStart)/CLOCKS_PER_SEC << endl;
  //MPI_Barrier( MPI_COMM_WORLD );
  #if RUNTIME == RUNTIME_VAL
  delete[] vec_calc_computation_time_X;
  delete[] vec_calc_computation_time_Y;
  delete[] vec_calc_thread_assign_X;
  delete[] vec_calc_thread_assign_Y;
  delete[] vec_calc_proc_assign_X;
  delete[] vec_calc_proc_assign_Y;
  delete[] chunk_gdvs_computation_time_X;
  delete[] chunk_gdvs_computation_time_Y;
  delete[] chunk_gdvs_proc_assign_X;
  delete[] chunk_gdvs_proc_assign_Y;
  #endif
  }
  Kokkos::finalize();
  MPI_Finalize(); 
  return 0;
}


KOKKOS_INLINE_FUNCTION void 
kokkos_Similarity_Metric_calculation_for_two_graphs(const matrix_type& graph1, 
                                                    const matrix_type& graph2, 
                                                    const Orbits& orbits, 
                                                    string graph_tag1, 
                                                    string graph_tag2,
						    int num_threads) {
  //clock_t out_tStart = clock();
  #if RUNTIME == RUNTIME_VAL
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
  kokkos_GDV_vector_calculation(graph1, graph1_GDV, orbits, "graph1", graph_counter, num_threads); 
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
if(rankm == 0)
  cout << "Done with GDV calculation for graph 1\n";
#endif

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
  kokkos_GDV_vector_calculation(graph2, graph2_GDV, orbits, "graph2", graph_counter, num_threads); 
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
if (rankm == 0)
  cout << "Done with GDV calculation for graph 2\n";
#endif

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

int m = 0;
int n = 0;
Kokkos::View<double**> sim_mat;

if (rankm == 0 && graph1_GDV.extent(0) < SIM_MAT_CUTOFF && graph2_GDV.extent(0) < SIM_MAT_CUTOFF) {
  m = graph1_GDV.extent(0);
  n = graph2_GDV.extent(0);
  //cout << "Creating similarity matrix with size " << m << "x" << n << endl;
  //Kokkos::View<double**> sim_mat("Similarity matrix", m, n);
  sim_mat = Kokkos::View<double**>("Similarity matrix", m, n);

  //cout << "Adding data to similarity matrix" << endl;
  Kokkos::MDRangePolicy<Kokkos::Rank<2>> mdrange({0,0}, {m,n});
  Kokkos::parallel_for("Compute sim_mat", mdrange, KOKKOS_LAMBDA(const int i, const int j) {
      auto g1_subview = Kokkos::subview(graph1_GDV, i, Kokkos::ALL);
      auto g2_subview = Kokkos::subview(graph2_GDV, j, Kokkos::ALL);
      sim_mat(i,j) = kokkos_GDV_distance_calculation(g1_subview, g2_subview);
  });
  //cout << "Constructed similarity matrix" << endl;
}
#ifdef DEBUG
if (rankm == 0) {
  cout << "Created similarity matrix\n";
}
#endif
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
  #if RUNTIME == RUNTIME_VAL
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
  #if RUNTIME == RUNTIME_VAL
  report_output_time = (double)(MPI_Wtime() - report_tStart);
  #endif
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
  unsigned long long sum = 0;
  for(int i=0; i<gdvm.extent(0); i++) {
    //sum += gdvm[i]*gdvm[i];
    sum += static_cast<double>(static_cast<int64_t>(gdvm[i])*static_cast<int64_t>(gdvm[i]));
  }
  return sqrt(sum);
}


void 
kokkos_GDV_vector_calculation(const matrix_type& graph, 
                              GDVs& graph_GDV, 
                              const Orbits& orbits, 
                              const char* graph_name, 
                              int graph_counter,
			      int num_threads) {
  // Set up parallelization               
  int comm_size, rankn;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankn);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  int tag = 11;
  int graph_size = graph.numRows();
//cout << "Starting GDV vector calculation\n";
//cout << "# of processes : " << comm_size << endl;
//  int graph_size = graph.size();
//  graph_GDV->clear();
//
//  for(int i=0; i<graph.size(); i++) {
//    vector<int> gdv_vec(GDV_LENGTH, 0);
//    graph_GDV->push_back(GDVMetric(graph[i].Row, gdv_vec));
//  }
  //graph_counter += 1;
  //cout << "Graph counter on rank " << rankn << " = " << graph_counter << endl;

#if RUNTIME == RUNTIME_VAL
  double process_ends_communication;
  double vec_calc_computation_start;
  double vec_calc_computation_end;
  double vec_calc_start = MPI_Wtime();
#endif

  Kokkos::View<unsigned long int*> num_combinations("Number of combinations", graph.numRows());
  //Kokkos::View<long long int*> num_combinations("Number of combinations", graph.numRows());
  //cout << "Initial Number of Combinations on Graph Size=" << graph.numRows() << " is : " << num_combinations(graph.numRows()-1) << endl;
  int k_interval = 1000000;

  Kokkos::View<int*> neighbor_scratch("Neighbors", graph.numRows());
  Kokkos::parallel_for(graph.numRows(), KOKKOS_LAMBDA(const int node) {
    int num_neighbors = 0;
    num_neighbors = EssensKokkos::get_num_neighbors(graph, node, 4, neighbor_scratch);
    for(int i=1; i<5; i++) {
      num_combinations(node) += get_num_combinations(num_neighbors, i);
      //cout << "Changing Number of Combinations at node " << node << " to: " << num_combinations(node) << " based on iteration " << i << " of get_num_neighbors=" << get_num_combinations(num_neighbors, i) << endl;
    }
  });
  //printf("Computed # of combinations\n");
  Kokkos::parallel_scan(num_combinations.extent(0), KOKKOS_LAMBDA(const int i, unsigned long int& update, const bool final) {
    const int val_i = num_combinations(i);
    if(final) {
      num_combinations(i) = update;
      //cout << "Setting Number of Combinations at node " << i << " to: " << num_combinations(i) << " based on update=" << update << endl;
    }
    update += val_i;
  });
  //cout << "Number of Combinations: " << num_combinations(graph.numRows()-1) << " on rank: " << rankn << endl;
  Kokkos::View<unsigned long int*> starts("Start indices", (num_combinations(graph.numRows()-1)/k_interval)+1);
  Kokkos::View<unsigned long int*> ends("Start indices", (num_combinations(graph.numRows()-1)/k_interval)+1);
  //Kokkos::View<long long int*> starts("Start indices", (num_combinations(graph.numRows()-1)/k_interval)+1);
  //Kokkos::View<long long int*> ends("Start indices", (num_combinations(graph.numRows()-1)/k_interval)+1);
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
  #ifdef DEBUG
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
  #endif

  if(comm_size == 1)
  {
//    Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>> policy(1, Kokkos::AUTO());
    dynamic_serial_work(graph, graph_GDV, orbits, starts, ends, num_threads, graph_counter);

#if RUNTIME == RUNTIME_VAL
    process_ends_communication = MPI_Wtime();
#endif
  }
  else if (rankn == 0) 
  {
    dynamic_master_work(graph_GDV, graph_size, comm_size);

#if RUNTIME == RUNTIME_VAL
      process_ends_communication = MPI_Wtime();
      //vec_calc_prior_gather = MPI_Wtime() - vec_calc_start + vec_calc_prior_gather;
#endif

      // Sort return vector
//      sort(graph_GDV->begin(), graph_GDV->end(), node_id_order());
//      #ifdef DEBUG
//        cout << "Constructed return GDV array" << endl;
//        for (i = 0; i < graph_GDV.size(); i++) {
//      	  cout << graph_GDV->at(i).node << ": ";
//      	  for (int j = 0; j < graph_GDV->at(i).GDV.size(); j++) {
//      	    cout << graph_GDV->at(i).GDV[j] << ", ";
//      	  }
//      	  cout << endl;
//      	}	
//      #endif
//    }
  }
  else // Instructions for work processes
  { 
    dynamic_worker_work(graph, orbits, num_threads, graph_counter);
//    graph_GDV->clear();
    #if RUNTIME == RUNTIME_VAL
        process_ends_communication = MPI_Wtime();
    #endif
  }

  //vec_calc_post_gather = MPI_Wtime() - vec_calc_start + vec_calc_post_gather;
#if RUNTIME == RUNTIME_VAL
  vec_calc_communication_time[graph_counter - 1] = process_ends_communication - vec_calc_start;
#endif
//  cout << "Communication time on process " << rankn << " for graph " << graph_counter << " = " << vec_calc_communication_time[graph_counter - 1] << endl;
  //vec_calc_computation_time = vec_calc_computation_end - vec_calc_computation_start + vec_calc_computation_time;

  #ifdef DEBUG
    cout << "Finished GDV Vector Calc on Rank: " << rankn << endl;
  #endif

}

template<class NeighborView, class IntView, class GraphView, class BoolView, class CounterView
>
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
                      CounterView& counter,
                      Kokkos::Experimental::ScatterView<int**> gdvMetrics_sa
                    )
{
  auto gdvMetrics = gdvMetrics_sa.access();
//  int combination_count = 0;
  int k_interval = 1000000;
  int num_neighbors = EssensKokkos::find_neighbours(node, graph, 4, neighbor_buff);
  auto neighbors = Kokkos::subview(neighbor_buff, std::pair<int,int>(0, num_neighbors));
#ifdef DEBUG
  printf("Allocated and found neighbors of %d: %d total neighbors\n", node, neighbors.size());
#endif
  int iter_counter = counter(team_member.team_rank());
  for (int node_count = 1; node_count < 5; node_count++)
  {
//printf("Subgraphs of size %d\n", node_count+1);
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
      int num_edges = EssensKokkos::kokkos_induced_subgraph(graph, combination, induced_sgraph);
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
      generator.kokkos_next(indices);

//      combination_count++;
      iter_counter += 1;
    }
//    printf("Number of combinations: %d for %d neighbors and %d node subgraphs\n", generator.get_num_comb(), num_neighbors, node_count+1);
  }
  counter(team_member.team_rank()) = iter_counter;
//  printf("Thread %d at iteration %d\n", team_member.team_rank(), iter_counter);
}


void dynamic_serial_work(const matrix_type& graph, GDVs& graph_GDV, const Orbits& orbits, Kokkos::View<unsigned long int*> starts, Kokkos::View<unsigned long int*> ends, int num_threads, int graph_counter) {
    int tag = 11;

    Kokkos::TeamPolicy<> policy(1, num_threads);
    using member_type = Kokkos::TeamPolicy<>::member_type;

    #ifdef DEBUG
      cout << "Team size: " << policy.team_size() << endl;
    #endif
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

#if RUNTIME == RUNTIME_VAL
    // Start measuring time for shared memory parallelism
    double chunk_gdvs_computation_start = MPI_Wtime();
#endif

    int i=0;
    for(i; i<starts.extent(0); i++) {
      Kokkos::parallel_for("Calcualte GDV bundle", policy, KOKKOS_LAMBDA(member_type team_member) {
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
#ifdef DEBUG
chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
#endif

#if RUNTIME == RUNTIME_VAL
        double vec_calc_computation_start = MPI_Wtime();
#endif
     
          kokkos_calculate_GDV(team_member, node, graph, orbits, neighbor_subview, indices_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, combination_counter,
                              metrics_sa);

#if RUNTIME == RUNTIME_VAL
        if (graph_counter == 1) {
          vec_calc_computation_time_X[node] = MPI_Wtime() - vec_calc_computation_start;
          vec_calc_thread_assign_X[node] = team_member.team_rank();
	  vec_calc_proc_assign_X[node] = 0;
        }
        if (graph_counter == 2) {
          vec_calc_computation_time_Y[node] = MPI_Wtime() - vec_calc_computation_start;
          vec_calc_thread_assign_Y[node] = team_member.team_rank();
	  vec_calc_proc_assign_Y[node] = 0;
        }
#endif

#ifdef DEBUG
          chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
          chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2-t1);
          printf("Done with node: %d, time: %f\n", node, time_span.count());
#endif
        });
      });

      Kokkos::Experimental::contribute(graph_GDV, metrics_sa);
      metrics_sa.reset();
    }

#if RUNTIME == RUNTIME_VAL
    // End time recording for shared memory parallelism
    if (graph_counter == 1) {
      chunk_gdvs_computation_time_X[0] = MPI_Wtime() - chunk_gdvs_computation_start;
      chunk_gdvs_proc_assign_X[0] = 0;
    }
    if (graph_counter == 2) {
      chunk_gdvs_computation_time_Y[0] = MPI_Wtime() - chunk_gdvs_computation_start;
      chunk_gdvs_proc_assign_Y[0] = 0;
    }
#endif

}

void dynamic_master_work(GDVs& graph_GDV, int graph_size, int comm_size) {
    int tag = 11;
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
      int chunk_index = 0;
      int node_chunk_package[2];
      int rcv_node;  // Corresponds to name of graph node recieved from worker rank
      for (i = 1; i < comm_size; i++) {
              send_node = (i-1) * CHUNK_SIZE;
	      node_chunk_package[0] = send_node;
	      node_chunk_package[1] = chunk_index;
        #ifdef DEBUG
          cout << "Sending node " << send_node << " to rank " << i << endl;
        #endif
        MPI_Send(node_chunk_package, 2, MPI_INT, i, tag, MPI_COMM_WORLD);
	chunk_index += 1;
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
//                      cout << "Sending node " << graph[send_node].Row << " to rank " << i << endl;
                cout << "Sending node " << send_node << " to rank " << i << endl;
              #endif
              node_chunk_package[0] = send_node;
              node_chunk_package[1] = chunk_index;
              MPI_Send(node_chunk_package, 2, MPI_INT, i, tag, MPI_COMM_WORLD);
	      chunk_index += 1;
              send_node += CHUNK_SIZE;
            } else { // Send termination
              flag = -1;
	      node_chunk_package[0] = flag;
              node_chunk_package[1] = chunk_index;
              MPI_Send(node_chunk_package, 2, MPI_INT, i, tag, MPI_COMM_WORLD);
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
//                free(gdv_array);
        }
      } while (rec_count < graph_size);
    }
}

void dynamic_worker_work(const matrix_type& graph, const Orbits& orbits, int num_threads, int graph_counter) {

    int comm_size, rankn;
    MPI_Comm_rank(MPI_COMM_WORLD, &rankn);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    int tag = 11;
    int node_name;
    int chunk_index;
    int node_chunk_array[2];
    MPI_Status worker_status;

    do {

      //cout << "About to do MPI_Recv on rank " << rankn << endl;
//cout << "Starting worker thread\n";
      MPI_Recv(node_chunk_array, 2, MPI_INT, 0, tag, MPI_COMM_WORLD, &worker_status);
      node_name = node_chunk_array[0];
      chunk_index = node_chunk_array[1];
//cout << "Got new batch of nodes starting at " << node_name << "\n";
      #ifdef DEBUG
        cout << "Recieved node " << node_name << " at rank " << rankn << endl;
      #endif
      if (node_name == -1) {
//#if RUNTIME == RUNTIME_VAL
//        process_ends_communication = MPI_Wtime();
//#endif
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
      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(1, num_threads);
//      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy((end-node_name)/4, 4);
//      policy.set_scratch_size(0, Kokkos::PerThread(128));
      using member_type = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;

#ifdef DEBUG
cout << "Team size: " << policy.team_size() << endl;
#endif
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
#ifdef DEBUG
chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
#endif

#if RUNTIME == RUNTIME_VAL
  double chunk_gdvs_computation_start = MPI_Wtime();
#endif

//      Kokkos::View<int[1]> node_counter("Node counter");
//      node_counter(0) = node_name;
      Kokkos::parallel_for("Calculate GDV", policy, KOKKOS_LAMBDA(member_type team_member) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, end-node_name), [=] (int n) {
//          int node = node_name + team_member.league_rank()*team_member.team_size() + team_member.team_rank();
//          int node = Kokkos::atomic_fetch_add(&node_counter(0), 1);
          int node = node_name + n;
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

#if RUNTIME == RUNTIME_VAL
        double vec_calc_computation_start = MPI_Wtime();
#endif

	kokkos_calculate_GDV(team_member, node, graph, orbits, neighbor_subview, indices_subview, combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, combination_counter,
                                metrics_sa);

#if RUNTIME == RUNTIME_VAL
        if (graph_counter == 1) {
          vec_calc_computation_time_X[node] = MPI_Wtime() - vec_calc_computation_start;
	  vec_calc_thread_assign_X[node] = rankn;
          vec_calc_proc_assign_X[node] = rankn;
        }
        if (graph_counter == 2) {
          vec_calc_computation_time_Y[node] = MPI_Wtime() - vec_calc_computation_start;
	  vec_calc_thread_assign_Y[node] = rankn;
          vec_calc_proc_assign_Y[node] = rankn;
        }
#endif

        });
      });

#ifdef DEBUG
chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2-t1);
printf("Done with node chunk: %d, time: %f\n", node_name, time_span.count());
std::cout << "Rank " << rankn << " finished chunk " << node_name << std::endl;
#endif

      Kokkos::Experimental::contribute(metrics, metrics_sa);
#if RUNTIME == RUNTIME_VAL
      // End time recording for shared memory parallelism
    if (graph_counter == 1) {
      chunk_gdvs_computation_time_X[chunk_index] = MPI_Wtime() - chunk_gdvs_computation_start;
      chunk_gdvs_proc_assign_X[chunk_index] = rankn;
    }
    if (graph_counter == 2) {
      chunk_gdvs_computation_time_Y[chunk_index] = MPI_Wtime() - chunk_gdvs_computation_start;
      chunk_gdvs_proc_assign_Y[chunk_index] = rankn;
    }
#endif

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
  int orbit_size = 2;
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
    for(int i=0; i<vector_line[1].size(); i++) {
      orbits.degree(orbit_counter, i) = vector_line[1][i];
    }
    for(int i=0; i<vector_line[2].size(); i++) {
      orbits.distance(orbit_counter, i) = vector_line[2][i];
    }
    orbit_counter++;
  }
  orbits.start_indices(6) = num_orbits;
}

#if RUNTIME == RUNTIME_VAL
void record_runtimes(const matrix_type& graphx, const matrix_type& graphy, double* time_buffer, int num_times, time_t now, int mpi_rank, int numtasks) {

  if (mpi_rank == 0) {

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
      myfile << month << "/" << day << "/" << year << "\t" << numtasks << "\t" << graph_name1 << "\t" << graph_name2 << "\t" << graphx.numRows() << "\t\t" << total_time_taken << " \n";
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
      myfile << i << " " << time_buffer[num_times*i+3] << " " << time_buffer[num_times*i+4] << " " << time_buffer[num_times*i+1] << " " << time_buffer[num_times*i+2] << " " << time_buffer[num_times*i] << " \n";
    }
    myfile.close();
  }

  string thread_computation_time_file_x = "runtime_data/runtimes_rec_for_nodes_x_" + to_string(mpi_rank) + ".txt";
  string thread_computation_time_file_y = "runtime_data/runtimes_rec_for_nodes_y_" + to_string(mpi_rank) + ".txt";
  ofstream myfile;
  myfile.open(thread_computation_time_file_x, ofstream::trunc);
  if (!myfile.is_open()) {
    cout << "INPUT ERROR:: Could not open the local time recording file\n";
  }
  if (myfile.is_open()) {
    for (int i = 0; i < graphx.numRows(); i++) {
      myfile << i << " " << vec_calc_proc_assign_X[i] << " " << vec_calc_thread_assign_X[i] << " " << vec_calc_computation_time_X[i] << endl;
    }
    myfile.close();
  }
  myfile.open(thread_computation_time_file_y, ofstream::trunc);
  if (!myfile.is_open()) {
    cout << "INPUT ERROR:: Could not open the local time recording file\n";
  }
  if (myfile.is_open()) {
    for (int i = 0; i < graphy.numRows(); i++) {
      myfile << i << " " << vec_calc_proc_assign_Y[i] << " " << vec_calc_thread_assign_Y[i] << " " << vec_calc_computation_time_Y[i] << endl;
    }
    myfile.close();
  }

  string chunk_computation_time_file_x = "runtime_data/runtimes_rec_for_chunks_x_" + to_string(mpi_rank) + ".txt";
  string chunk_computation_time_file_y = "runtime_data/runtimes_rec_for_chunks_y_" + to_string(mpi_rank) + ".txt";
  myfile.open(chunk_computation_time_file_x, ofstream::trunc);
  if (!myfile.is_open()) {
    cout << "INPUT ERROR:: Could not open the local time recording file\n";
  }
  if (myfile.is_open()) {
    for (int i = 0; i < ((graphx.numRows()/CHUNK_SIZE)+1); i++) {
      myfile << i << " " << chunk_gdvs_proc_assign_X[i] << " " << chunk_gdvs_computation_time_X << " " << endl;
    }
    myfile.close();
  }
  myfile.open(chunk_computation_time_file_y, ofstream::trunc);
  if (!myfile.is_open()) {
    cout << "INPUT ERROR:: Could not open the local time recording file\n";
  }
  if (myfile.is_open()) {
    for (int i = 0; i < ((graphy.numRows()/CHUNK_SIZE)+1); i++) {
      myfile << i << " " << chunk_gdvs_proc_assign_Y[i] << " " << chunk_gdvs_computation_time_Y << " " << endl;
    }
    myfile.close();
  }

}
#endif

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

/*
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
    if ( (m < SIM_MAT_CUTOFF) || (n < SIM_MAT_CUTOFF) ) {
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
//      	  #ifdef DEBUG
//            cout << "Recieved GDV for node " << rcv_node << " from rank " << i << ": " << endl;
//            cout << gdv_array[GDV_LENGTH] << ": ";
//            for (int j = 0; j < GDV_LENGTH; j++) {
//              cout << gdv_array[j] << ", ";
//            }
//            cout << endl;
//          #endif
      
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
#endif

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
#ifdef RUNTIME
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
  vector<int> neighbours;
  gdvf.find_neighbours(node,Graph,4,&neighbours);
  int set[neighbours.size()]; 
  std::copy( neighbours.begin(), neighbours.end(), set );
  int numElements = sizeof(set)/sizeof(set[0]);
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
*/


