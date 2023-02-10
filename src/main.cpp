//#define AUTO_CHECKPOINT
#include "deduplicator.hpp"
#include "class_definitions.hpp"
#include "combinations.hpp"
#include "kokkos_functions.hpp"
#include "dirty_page_tracking.hpp"
#include "hash_change_detection.hpp"
#include "update_pattern_analysis.hpp"
#include "io.hpp"
//#include "hash_functions.hpp"
//#include "deduplication.hpp"
#include "utils.hpp"
#include <time.h>
#include <stdlib.h>
#include <ctime>
#include <mpi.h>
#include <fstream>
#include <chrono>
#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include <Kokkos_ScatterView.hpp>
#include <Kokkos_UnorderedMap.hpp>

#define GDV_LENGTH 73
#define NUM_THREADS Kokkos::AUTO
//#define NUM_THREADS 1
#define MAX_SUBGRAPH 5
//#define ANALYSIS
//#define HASH_ANALYSIS
//#define SCAN_ANALYSIS
//#define DEBUG
//#define DIRTY_PAGE_TRACKING
//#define AUTO_CHECKPOINT
//#define OUTPUT_MATRIX
//#define OUTPUT_GDV
//#define LAYOUT Kokkos::LayoutLeft
#define DETAILED_TIMERS

//#define DEDUPLICATION
//#define DEDUP_CASE1

#ifdef AUTO_CHECKPOINT
#include <resilience/Resilience.hpp>
#include <resilience/CheckpointFilter.hpp>
#include <veloc.h>
#endif

void compute_gdvs(const matrix_type& graph, const Orbits& orbits, GDVs& graph_GDV, string name);
void kokkos_similarity_metric_calculation_for_two_graphs(const matrix_type&, const matrix_type&, const Orbits&, string, string);
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
//template<class NeighborView, class IntView, class GraphView, class BoolView, class Scatter>
template<class IntView, class GraphView, class BoolView, class Scatter>
KOKKOS_INLINE_FUNCTION void 
kokkos_calculate_GDV2(Kokkos::TeamPolicy<>::member_type team_member,
                      int node, 
                      int node_count,
                      const matrix_type& graph, 
                      const Orbits& orbits, 
//                      const NeighborView& neighbors,
                      const int n_neighbors,
                      IntView& combination_view,
                      IntView& sgraph_distance_signature,
                      IntView& sgraph_degree_signature,
                      GraphView& induced_subgraph,
                      BoolView& visited,
                      IntView& queue,
                      IntView& distance,
                      Scatter gdvMetrics_sa
                    );

// Define variables for keeping track of time for load imbalancing tests.
double total_time_taken;
double vec_calc_communication_time[2] = {};
double pre_process_time;
double report_output_time;
int vec_calc_avg_node_deg;
uint64_t CHECKPOINT_INTERVAL;
int CHUNK_SIZE;
uint64_t MAX_INTERVALS;
uint64_t TIME_INTERVAL;
//uint64_t NUM_CHECKPOINTS=4;

using namespace std;

DedupMode dedup_mode;

int main(int argc, char *argv[]) {
  int numtasks, rank;//, dest, source, rc, count, tag=0;

  MPI_Init(&argc,&argv);
  Kokkos::initialize(argc, argv);
  {

  time_t now = time(0);
  double pre_tStart = MPI_Wtime();
  double start_time, end_time;
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

  ifstream the_file3 ( argv[4] );
  if (!the_file3.is_open() ) {
    cout << "INPUT ERROR:: Could not open the time recording file\n";
  }

  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) 
    printf("Comm size: %d\n", numtasks);

  if(rank == 0) {
    start_time = MPI_Wtime();
  }

  CHECKPOINT_INTERVAL = strtoull(argv[5], NULL, 0);
  sscanf(argv[6], "%d", &CHUNK_SIZE);

  sscanf(argv[7], "%lu", &MAX_INTERVALS);

  sscanf(argv[8], "%lu", &TIME_INTERVAL);


  dedup_mode = get_mode(argc, argv);

  Orbits k_orbits;
  kokkos_readin_orbits(&the_file2, k_orbits);

  matrix_type graphx, graphy;
  readin_graph(&the_file0, graphx);
//  readin_graph(&the_file1, graphy);

  if(rank == 0) {
    end_time = MPI_Wtime();
    printf("Time spent reading data: %lf\n", end_time-start_time);
  }

//  MPI_Status Stat;   // required variable for receive routines

  // Get data on graph names
  string delimiter = "/";
  string graph_name1 = argv[1];
  string graph_name2 = argv[2];
  
  uint64_t pos1 = 0;
  string token1;
  while ((pos1 = graph_name1.find(delimiter)) != string::npos) {
    token1 = graph_name1.substr(0, pos1);
    //cout << token << endl;
    graph_name1.erase(0, pos1 + delimiter.length());
  }
  
  uint64_t pos2 = 0;
  string token2;
  while ((pos2 = graph_name2.find(delimiter)) != string::npos) {
    token2 = graph_name2.substr(0, pos2);
    graph_name2.erase(0, pos2 + delimiter.length());
  }

  #ifdef DEBUG
    cout << "Starting Similarity Metric Calculation on Rank: " << rank << endl;
  #endif

  if(rank == 0) {
    cout << "Graph 1: # of vertices: " << graphx.numRows() << endl;
//    cout << "Graph 2: # of vertices: " << graphy.numRows() << endl;
  }

  pre_process_time = (double)(MPI_Wtime() - pre_tStart);
  if(rank == 0) {
    start_time = MPI_Wtime();
  }

  // Perform Similarity Calculations
//  kokkos_similarity_metric_calculation_for_two_graphs(graphx, graphy, k_orbits, graph_name1, graph_name2);

  using namespace std::chrono;
  high_resolution_clock::time_point start_time0 = high_resolution_clock::now();
  int num_orbits = k_orbits.degree.extent(0);
  GDVs gdvs("GDVs", graphx.numRows(), num_orbits);
  GDVs::HostMirror graph1_gdvs = Kokkos::create_mirror_view(gdvs);
  compute_gdvs(graphx, k_orbits, gdvs, graph_name1);
  Kokkos::deep_copy(graph1_gdvs, gdvs);
  high_resolution_clock::time_point end_time0 = high_resolution_clock::now();
  printf("Time spent on graph 1: %f\n", duration_cast<duration<double>>(end_time0 - start_time0));

//  readin_graph(&the_file1, graphy);
//  if(rank == 0) {
////    cout << "Graph 1: # of vertices: " << graphx.numRows() << endl;
//    cout << "Graph 2: # of vertices: " << graphy.numRows() << endl;
//  }
//  CHECKPOINT_INTERVAL = strtoull(argv[5], NULL, 0);
//  sscanf(argv[6], "%d", &CHUNK_SIZE);
//  dedup_mode = get_mode(argc, argv);
//  start_time0 = high_resolution_clock::now();
//  Kokkos::resize(gdvs, graphy.numRows(), num_orbits);
//  Kokkos::deep_copy(gdvs, 0);
//  compute_gdvs(graphy, k_orbits, gdvs, graph_name2);
//  GDVs::HostMirror graph2_gdvs = Kokkos::create_mirror_view(gdvs);
//  Kokkos::deep_copy(graph2_gdvs, gdvs);
//  end_time0 = high_resolution_clock::now();
//  printf("Time spent on graph 2: %f\n", duration_cast<duration<double>>(end_time0 - start_time0));

#ifdef OUTPUT_GDV    
  write_gdvs(graph1_gdvs, graph_name1, graph2_gdvs, graph_name2);
#endif

  if(rank == 0) {
    end_time = MPI_Wtime();
    printf("Rank %d spent %lf s on calculating GDVs\n", rank, end_time-start_time);
  }

  #ifdef DEBUG
    cout << "Finished Similarity Metric Calculation on Rank: " << rank << endl;
  #endif

  // Set up for and Perform Runtime Management and Gather
  double* time_buff = NULL;
  int num_times = 5;
  double send_times[num_times];
  send_times[0] = total_time_taken;
  send_times[1] = vec_calc_communication_time[0];
  send_times[2] = vec_calc_communication_time[1];
  send_times[3] = pre_process_time;
  send_times[4] = report_output_time;
  if (rank == 0) {
    time_buff = (double *)malloc(numtasks*num_times*sizeof(double));
  }
  MPI_Gather(send_times, num_times, MPI_DOUBLE, time_buff, num_times, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  #ifdef DEBUG
    cout << "Finished time gathering and prepped for time-keeping on Rank: " << rank << endl;
  #endif

  start_time = MPI_Wtime();
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
      uint64_t team_size = team_policy.team_size();
      if(team_size > 0) {
        myfile << month << "/" << day << "/" << year << "\t" << numtasks << "\t" << team_policy.team_size() << "\t" << graph_name1 << "\t" << graph_name2 << "\t" << graphx.numRows() << "\t\t" << total_time_taken << " \n";
      } else {
        myfile << month << "/" << day << "/" << year << "\t" << numtasks << "\t" << "AUTO" << "\t" << graph_name1 << "\t" << graph_name2 << "\t" << graphx.numRows() << "\t\t" << total_time_taken << " \n";
      }
      myfile.close();
    } else { 
      cout << "Out File Did Not Open";
    }

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

  end_time = MPI_Wtime();
  if(rank == 0) {
    printf("Time spent on output: %lfs\n", end_time-start_time);
  }

  if(rank == 0)
    printf("Time spent overall: %lf\n", MPI_Wtime()-pre_tStart);
  }
  Kokkos::finalize();
  MPI_Finalize(); 
  return 0;
}

void compute_gdvs(const matrix_type& graph, const Orbits& orbits, GDVs& graph_GDV, string name) {
  MPI_Barrier( MPI_COMM_WORLD );
  int rankm;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankm);
  kokkos_GDV_vector_calculation(graph, graph_GDV, orbits, name.c_str(), 0); 
  Kokkos::fence();
  MPI_Barrier(MPI_COMM_WORLD);
//#ifdef DEBUG
//  GDVs::HostMirror graph_GDV_host = Kokkos::create_mirror_view(graph_GDV);
//  Kokkos::deep_copy(graph_GDV_host, graph_GDV);
//  Kokkos::fence();
//  if(rankm == 0) {
//    cout << "Post kokkos_GDV_vector_calculation: graph " << name << endl;
//    for(int i=0; i<graph_GDV.extent(0); i++) {
//      cout << "Node " << i << ": ";
//      for(int j=0; j<graph_GDV_host.extent(1); j++) {
//        cout << graph_GDV_host(i,j) << " ";
//      }
//      cout << endl;
//    }
//    cout << endl;
//  }
//#endif
}

void 
kokkos_similarity_metric_calculation_for_two_graphs(const matrix_type& graph1, 
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
  GDVs::HostMirror graph1_GDV_host = Kokkos::create_mirror_view(graph1_GDV);
  GDVs::HostMirror graph2_GDV_host = Kokkos::create_mirror_view(graph1_GDV);

  int rankm, numtasksm;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankm);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasksm);

  double start;
  if(rankm == 0)
    start = MPI_Wtime();
  int graph_counter = 1;
  kokkos_GDV_vector_calculation(graph1, graph1_GDV, orbits, graph_tag1.c_str(), graph_counter); 
  Kokkos::fence();
  MPI_Barrier(MPI_COMM_WORLD);
  if(rankm == 0) {
    printf("Time spent calculating GDVs of graph 1: %lf\n", MPI_Wtime()-start);
    cout << "Done with GDV calculation for graph 1\n";
  }

//#ifdef DEBUG
//  Kokkos::deep_copy(graph1_GDV_host, graph1_GDV);
//  Kokkos::fence();
//  if(rankm == 0) {
//    cout << "Post kokkos_GDV_vector_calculation: graph 1" << endl;
//    for(int i=0; i<graph1_GDV.extent(0); i++) {
//      cout << "Node " << i << ": ";
//      for(int j=0; j<graph1_GDV_host.extent(1); j++) {
//        cout << graph1_GDV_host(i,j) << " ";
//      }
//      cout << endl;
//    }
//    cout << endl;
//  }
//#endif

//#ifdef AUTO_CHECKPOINT
//  VELOC_Checkpoint_wait();
//#endif

  if(rankm == 0)
    start = MPI_Wtime();
  graph_counter = 2;
  kokkos_GDV_vector_calculation(graph2, graph2_GDV, orbits, graph_tag2.c_str(), graph_counter); 
  MPI_Barrier(MPI_COMM_WORLD);
  if(rankm == 0) {
    printf("Time spent calculating GDVs of graph 1: %lf\n", MPI_Wtime()-start);
    cout << "Done with GDV calculation for graph 2\n";
  }

//#ifdef DEBUG
//  Kokkos::deep_copy(graph2_GDV_host, graph2_GDV);
//  if(rankm == 0) {
//    cout << "Post kokkos_GDV_vector_calculation: graph 2" << endl;
//    for(int i=0; i<graph2_GDV.extent(0); i++) {
//      cout << "Node " << i << ": ";
//      for(int j=0; j<graph2_GDV.extent(1); j++) {
//        cout << graph2_GDV_host(i,j) << " ";
//      }
//      cout << endl;
//    }
//    cout << endl;
//  }
//#endif

#ifdef OUTPUT_MATRIX
  Kokkos::fence();
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
  Kokkos::fence();
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
    Kokkos::View<double**>::HostMirror sim_mat_host = Kokkos::create_mirror_view(sim_mat);
    Kokkos::deep_copy(sim_mat_host, sim_mat);
    Kokkos::fence();
    string filename = "out_similarity_matrix.txt";
    myfile.open(filename);
    for(int i=1; i<m;i++) {
      for(int j=1;j<n;j++) {
        myfile<<" { "<<sim_mat_host(i,j)<<" } ";
      }
      myfile<<"||"<<endl;
    }
    myfile.close();
#endif
    
#ifdef OUTPUT_GDV    
    // Print out GDVs into files
    Kokkos::deep_copy(graph1_GDV_host, graph1_GDV);
    Kokkos::deep_copy(graph2_GDV_host, graph2_GDV);
    string gdv_file1 = "out_gdv_1_" + graph_tag1 + ".txt";
    string gdv_file2 = "out_gdv_2_" + graph_tag2 + ".txt";
    myfile.open(gdv_file1, ofstream::trunc);
    for (int i = 0; i < graph1_GDV_host.extent(0); i++) {
      for (int j = 0; j< graph1_GDV_host.extent(1); j++) {
        myfile << graph1_GDV_host(i,j) << " ";
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
  for(uint64_t i=0; i<gdvm.extent(0); i++) {
    sum += static_cast<double>(static_cast<int64_t>(gdvm[i])*static_cast<int64_t>(gdvm[i]));
  }
  return sqrt(sum);
}

void 
kokkos_GDV_vector_calculation(const matrix_type& graph, 
                              GDVs& graph_GDV, 
                              const Orbits& orbits, 
                              const char* graph_name, 
                              int graph_counter) {
  // Set up MPI parallelization               
  int comm_size, rankn;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankn);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
//  int tag = 11;
  uint32_t graph_size = graph.numRows();
  std::string label(graph_name);
  label += "_Rank";
  label += std::to_string(rankn);

#ifdef DEBUG
  cout << "Starting GDV vector calculation\n";
  cout << "# of processes : " << comm_size << endl;
#endif

  double process_ends_communication;
  double vec_calc_start = MPI_Wtime();

  Kokkos::View<uint64_t*> num_combinations("Number of combinations", graph.numRows());
  Kokkos::View<uint32_t*> num_neighbors("Number of neighbors", graph.numRows());
  Kokkos::View<uint32_t*>::HostMirror num_neighbors_h = Kokkos::create_mirror_view(num_neighbors);
  uint64_t k_interval = CHECKPOINT_INTERVAL;

  typedef Kokkos::DefaultExecutionSpace::scratch_memory_space ScratchSpace;
  typedef Kokkos::View<int*, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_int_view;
  typedef Kokkos::View<int**, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_subgraph_view;
  typedef Kokkos::View<bool*, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_bool_view;
#ifdef DEBUG
  printf("Starting combination counter kernel\n");
#endif

  uint32_t vertices_per_proc = graph.numRows()/comm_size;
  uint64_t start_offset = rankn*vertices_per_proc;
  uint64_t end_offset = start_offset+vertices_per_proc;
  if(end_offset > graph.numRows())
    end_offset = graph.numRows();
  int* recv_count = new int[comm_size];
  int* displs = new int[comm_size];
  int running_count = 0;
  for(uint64_t i=0; i<comm_size-1; i++) {
    recv_count[i] = vertices_per_proc;
    displs[i] = running_count;
    running_count += vertices_per_proc;
  }
  recv_count[uint64_t (comm_size-1)] = graph.numRows()-(comm_size-1)*vertices_per_proc;
  displs[comm_size-1] = running_count;
#ifdef DEBUG
  for(int i=0; i<comm_size; i++) {
    printf("(%d,%d) ", recv_count[i], displs[i]);
  }
  printf("\n");
  printf("(%u,%u)\n", start_offset, end_offset);
#endif

  chrono::steady_clock::time_point count_comb_beg = chrono::steady_clock::now();

//  uint64_t comb_count_scratch = sizeof(scratch_int_view) + sizeof(int)*graph.numRows();
//  Kokkos::TeamPolicy<> comb_count_policy = Kokkos::TeamPolicy<>(1, NUM_THREADS).set_scratch_size(1, Kokkos::PerThread(comb_count_scratch));
//  using team_member_type = Kokkos::TeamPolicy<>::member_type;
//  Kokkos::parallel_for("Get # of combinations", comb_count_policy, KOKKOS_LAMBDA(team_member_type team_member) {
//    scratch_int_view neighbor_subview(team_member.thread_scratch(1), graph.numRows());
//    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, end_offset-start_offset), [=] (int n) {
//      uint64_t node = start_offset+n;
//      uint32_t n_neighbors = 0;
//      n_neighbors = EssensKokkos::get_num_neighbors(graph, node, 4, neighbor_subview);
//      num_neighbors(node) = n_neighbors;
//      for(int i=1; i<5; i++) {
//        num_combinations(node) += get_num_combinations(n_neighbors, i);
//      }
//    });
//  });

  uint32_t max_degree = 0;
  Kokkos::parallel_reduce("Find degree", Kokkos::RangePolicy<>(0, graph.numRows()), KOKKOS_LAMBDA(const uint32_t node, uint32_t& sum) {
    uint32_t degree = graph.row(node).length;
    num_neighbors(node) = degree;
    if(degree > sum)
      sum = degree;
  }, Kokkos::Max<uint32_t>(max_degree));
  printf("Max degree: %u\n", max_degree);
  uint32_t max_neighborhood = max_degree + max_degree*max_degree + max_degree*max_degree*max_degree + max_degree*max_degree*max_degree*max_degree;
  if(max_neighborhood > graph.numRows())
    max_neighborhood = graph.numRows();
  printf("Approximate max degree: %u\n", max_neighborhood);

  uint32_t comb_n_leagues = 1;
  uint64_t comb_count_scratch = sizeof(scratch_int_view) + sizeof(int)*max_neighborhood;
  Kokkos::TeamPolicy<> comb_count_policy = Kokkos::TeamPolicy<>(comb_n_leagues, NUM_THREADS).set_scratch_size(1, Kokkos::PerThread(comb_count_scratch));
  using team_member_type = Kokkos::TeamPolicy<>::member_type;
  Kokkos::parallel_for("Get # of combinations", comb_count_policy, KOKKOS_LAMBDA(team_member_type team_member) {
    uint32_t league_rank = team_member.league_rank();
if(league_rank == 0 && team_member.team_rank() == 0) 
  printf("Counting combinations with %d threads\n", team_member.team_size());
    uint32_t per_league = (end_offset-start_offset)/comb_n_leagues;
    if(per_league*comb_n_leagues < (end_offset-start_offset))
      per_league += 1;
    if(start_offset+per_league*(league_rank+1) > end_offset)
      per_league = end_offset - (start_offset+per_league*league_rank);
    scratch_int_view neighbor_subview(team_member.thread_scratch(1), max_neighborhood);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, per_league), [=] (int n) {
      uint64_t node = start_offset+per_league*league_rank+n;
      uint32_t n_neighbors = 0;
      n_neighbors = EssensKokkos::get_num_neighbors(graph, node, 4, neighbor_subview);
      num_neighbors(node) = n_neighbors;
      for(int i=1; i<5; i++) {
        num_combinations(node) += get_num_combinations(n_neighbors, i);
      }
    });
  });
  Kokkos::fence();
  Kokkos::deep_copy(num_neighbors_h, num_neighbors);
  chrono::steady_clock::time_point count_comb_end = chrono::steady_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(count_comb_end-count_comb_beg);
  printf("Time spent counting combinations: %f\n", time_span.count());

#ifdef DEBUG
  printf("Counted number of combinations for each vertex\n");
#endif

  // Find the maximum number of neighbors within distance 4
  uint32_t max_neighbors = 0;
  Kokkos::parallel_reduce("Find max neighbors", graph.numRows(), KOKKOS_LAMBDA(const uint32_t& x, uint32_t& deg) {
    if(num_neighbors(x) > deg)
      deg = num_neighbors(x);
  }, Kokkos::Max<uint32_t>(max_neighbors));
  uint32_t n_neighbors;
  MPI_Allreduce(&max_neighbors, &n_neighbors, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  max_neighbors = n_neighbors;
  if(rankn == 0) 
    printf("Max degree: %d\n", max_neighbors);

  // Get the total # of combinations of potential subgraphs
  uint64_t total_combinations = 0;
  uint64_t local_combinations = 0;
  Kokkos::parallel_reduce("Number of combinations", end_offset-start_offset, 
  KOKKOS_LAMBDA(const int i, uint64_t& update) {
    update += num_combinations(i+start_offset);
  }, local_combinations);
  uint64_t total_neighbors = 0;
  uint64_t local_neighbors = 0;
  Kokkos::parallel_reduce("Number of neighbors", end_offset-start_offset, KOKKOS_LAMBDA(const uint64_t i, uint64_t& update) {
    update += num_neighbors(i);
  }, local_neighbors);
  Kokkos::View<uint64_t*> send_comb("Send buffer", end_offset-start_offset);
  Kokkos::View<uint64_t*>::HostMirror send_comb_host = Kokkos::create_mirror_view(send_comb);
  auto s_view = Kokkos::subview(num_combinations, Kokkos::make_pair(start_offset, end_offset));
  Kokkos::deep_copy(send_comb, s_view);
  Kokkos::deep_copy(send_comb_host, send_comb);
  Kokkos::View<uint64_t*>::HostMirror num_combinations_host = Kokkos::create_mirror_view(num_combinations);
  Kokkos::deep_copy(num_combinations_host, num_combinations);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allgatherv(send_comb_host.data(), 
                 end_offset-start_offset, 
                 MPI_UNSIGNED_LONG, 
                 num_combinations_host.data(), 
                 recv_count, 
                 displs, 
                 MPI_UNSIGNED_LONG, 
                 MPI_COMM_WORLD);
  delete[] recv_count, displs;
  MPI_Allreduce(&local_combinations, &total_combinations, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_neighbors, &total_neighbors, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  Kokkos::deep_copy(num_combinations, num_combinations_host);
Kokkos::fence();

  // Split combinations into intervals
  uint64_t num_intervals = total_combinations/k_interval;
  if(num_intervals*k_interval < total_combinations) {
    num_intervals++;
  }
  if(rankn == 0) {
    printf("Computed # of combinations: (%lu) split into %lu groups\n", total_combinations, num_intervals);
  }
#ifdef DEBUG
  printf("Rank %d: ", rankn);
  for(uint32_t i=0; i<num_intervals+1; i++) {
//    printf("%ld ", starts_host(i));
    printf("%ld ", get_node_from_combinations(graph, num_combinations_host, k_interval, i));
  }
  printf("\n");
#endif

// Calculate GDVs
#ifdef DEBUG
    printf("Multiple ranks\n");
#endif
    Kokkos::deep_copy(graph_GDV, 0);
    uint64_t intervals_per_rank = num_intervals/comm_size;
    uint64_t remaining_intervals = num_intervals%comm_size;
    if(rankn < remaining_intervals)
      intervals_per_rank += 1;
    uint64_t start_index = 0;
    if(rankn < remaining_intervals) {
      start_index = intervals_per_rank*rankn;
    }
    if(rankn >= remaining_intervals) {
      start_index = ((intervals_per_rank+1)*remaining_intervals) + 
                    (intervals_per_rank*(rankn-remaining_intervals));
    }
    if(start_index+intervals_per_rank > num_intervals) {
      intervals_per_rank = num_intervals - start_index;
    }
    int offset = 0;
#ifdef DEBUG
printf("Rank %d: start_index: %llu, intervals_per_rank: %llu\n", rankn, start_index, intervals_per_rank);
#endif
//    Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(1, NUM_THREADS);
    using member_type = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;

    typedef Kokkos::DefaultExecutionSpace::scratch_memory_space ScratchSpace;
    typedef Kokkos::View<int*, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_int_view;
    typedef Kokkos::View<int**, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_subgraph_view;
    typedef Kokkos::View<bool*, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_bool_view;

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

    const int num_leagues = 1;

    Deduplicator<MD5Hash> deduplicator(CHUNK_SIZE);
    std::vector<Kokkos::View<uint8_t*>::HostMirror> diffs;

    double prior_time = MPI_Wtime();
    uint32_t chkpt_counter = 0;

    uint64_t end_interval = MAX_INTERVALS;
    if(MAX_INTERVALS==0)
      end_interval = intervals_per_rank;
    printf("Rank %d: Chunk size: %d, Chkpt Interval: %lu, Max Intervals: %lu, Time interval: %lu\n", rankn, CHUNK_SIZE, CHECKPOINT_INTERVAL, MAX_INTERVALS, TIME_INTERVAL);
    while(offset < intervals_per_rank && (chkpt_counter < end_interval)) {
#ifdef DETAILED_TIMERS
      chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
#endif
      uint32_t chunk_idx = start_index+offset;
#ifdef DEBUG
printf("Rank %d: Chunk idx: %u, Chunk combination start: %zu\n", rankn, chunk_idx, chunk_idx*k_interval);
#endif

#ifdef AUTO_CHECKPOINT
KokkosResilience::checkpoint(*ctx, label, offset, [=] () mutable {
#endif
      int prev_node = get_node_from_combinations(graph, num_combinations_host, k_interval, (chunk_idx-1));
      int curr_node = get_node_from_combinations(graph, num_combinations_host, k_interval, chunk_idx);
      int next_node = get_node_from_combinations(graph, num_combinations_host, k_interval, (chunk_idx+1));
#ifdef DEBUG
      printf("Rank %d: Prev node: %d, Curr node: %d, Next node: %d\n", rankn, prev_node, curr_node, next_node);
#endif
      if((chunk_idx+1)*k_interval > total_combinations) {
        next_node = graph.numRows();
      }
      bool multi_node = (chunk_idx+1 == num_intervals+1 || next_node != curr_node) && 
                        (chunk_idx==0 || prev_node != curr_node);
#ifdef DEBUG
      printf("Rank %d: Multi node chunk: %d\n", rankn, multi_node);
#endif
      if(multi_node) {
#ifdef DEBUG
printf("Multi node\n");
#endif
        int start_node = curr_node;
        int end_node;
        if(chunk_idx+1 == num_intervals+1) {
          end_node = graph.numRows();
        } else {
          end_node = next_node;
        }

        int per_league = (end_node-start_node)/num_leagues;
        if(per_league*num_leagues < end_node-start_node)
          per_league += 1;
#ifdef DEBUG
printf("Per league: %d\n", per_league);
#endif

        uint64_t scratch_bytes = 8*sizeof(scratch_int_view) + sizeof(scratch_subgraph_view) + 
                                 sizeof(scratch_bool_view) + sizeof(bool)*5 + 
                                 (2*max_neighbors+61)*sizeof(int);
#ifdef DEBUG
printf("%lu bytes\n", scratch_bytes);
#endif
        Kokkos::TeamPolicy<> team_policy = Kokkos::TeamPolicy<>(per_league, NUM_THREADS).set_scratch_size(1, Kokkos::PerThread(scratch_bytes));
#ifdef DETAILED_TIMERS
        chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
#endif
        Kokkos::parallel_for("Calculate GDV", team_policy, KOKKOS_LAMBDA(member_type team_member) {
          int32_t node = start_node+team_member.league_rank();
          scratch_int_view neighbor_subview(team_member.thread_scratch(1), max_neighbors);
          scratch_int_view indices_subview(team_member.thread_scratch(1), 5);
          scratch_int_view combination_subview(team_member.thread_scratch(1), 5);
          scratch_int_view sgraph_dist_subview(team_member.thread_scratch(1), 6);
          scratch_int_view sgraph_deg_subview(team_member.thread_scratch(1), 5);
          scratch_int_view queue_subview(team_member.thread_scratch(1), 5);
          scratch_int_view distance_subview(team_member.thread_scratch(1), 5);
          scratch_subgraph_view subgraph_subview(team_member.thread_scratch(1), 5,5);
          scratch_bool_view visited_subview(team_member.thread_scratch(1), 5);
          scratch_int_view scratch_view(team_member.thread_scratch(1), max_neighbors);
          int n_neighbors = EssensKokkos::find_neighbours(node, graph, 4, neighbor_subview);
          for(int node_count=1; node_count<5; node_count++) {
            uint64_t n_combinations = get_num_combinations(n_neighbors, node_count);
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_combinations), [=] (uint64_t idx) {
              uint64_t combination_num = idx;
              combination_from_position(scratch_view, combination_num, n_neighbors, node_count);
              for(int j=0; j<node_count; j++) {
                  combination_subview(j) = neighbor_subview(scratch_view(j));
              }
              combination_subview(node_count) = node;
              kokkos_calculate_GDV2(team_member, node, node_count, graph, orbits, num_neighbors(node), combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, graph_GDV);
            });
          }
        });

Kokkos::fence();
#ifdef DEBUG
printf("Rank %d done with chunk %d\n", rankn, chunk_idx);
printf("Exited branch0\n");
#endif
      } else {
#ifdef DEBUG
printf("Rank %d: Split node\n", rankn);
#endif
Kokkos::fence();
        Kokkos::View<int*> neighbor_scratch("Neighbor scratch", graph.numRows());
        Kokkos::View<int[1]> num_neighbors_temp_d("Num neighbors temp");
        
        Kokkos::parallel_for(1, KOKKOS_LAMBDA(const uint32_t _i) {
          num_neighbors_temp_d(0) = EssensKokkos::find_neighbours(curr_node, graph, 4, neighbor_scratch);
#ifdef DEBUG
          printf("Rank %d: Node %d has %d neighbors within a distance of 4\n", rankn, curr_node, num_neighbors_temp_d(0));
#endif
        });
        auto num_neighbors_temp_h = Kokkos::create_mirror_view(num_neighbors_temp_d);
        Kokkos::deep_copy(num_neighbors_temp_h, num_neighbors_temp_d);
        int n_neighbors = num_neighbors_temp_h(0);
Kokkos::fence();
#ifdef DEBUG
printf("Rank %d: Number of neighbors for %d: %d\n", rankn, curr_node, n_neighbors);
#endif

        uint64_t start_combination;
        uint64_t end_combination;
        uint64_t start_comb_subgraph[5] = {0,0,0,0,0};
        uint64_t end_comb_subgraph[5] = {0,0,0,0,0};
        start_comb_subgraph[0] = 0;
        end_comb_subgraph[0] = 0;
        int start_chunk = 0;
#ifdef DEBUG
printf("Rank: %d: Initialized split data\n", rankn);
#endif
        Kokkos::parallel_reduce("Calculate start chunk", Kokkos::RangePolicy<>(0, chunk_idx), KOKKOS_LAMBDA(const uint64_t j, int& update) {
          int u = get_node_from_combinations(graph, num_combinations, k_interval, j);
          if(u == curr_node)
            update += 1;
        }, start_chunk);
        Kokkos::fence();
#ifdef DEBUG
printf("Rank: %d: Start split section, start_chunk: %d\n", rankn, start_chunk);
#endif
        bool first_chunk=false, middle_chunk=false, last_chunk=false;
        first_chunk = ((chunk_idx == 0) || 
                        prev_node != curr_node) && 
                        next_node == curr_node;
        if(!first_chunk) {
          last_chunk = prev_node == curr_node && 
                       next_node != curr_node;
          middle_chunk = prev_node == curr_node && 
                         next_node == curr_node;
        }
#ifdef DEBUG
printf("Rank %d: Chunk (node %d) is first: %d, middle: %d, last: %d\n", rankn, curr_node, first_chunk, middle_chunk, last_chunk);
#endif
        if(last_chunk) {
          // Last chunk
          start_combination = start_chunk*k_interval;
          end_combination = num_combinations_host(curr_node);
#ifdef DEBUG
printf("Rank %d: Last chunk: start comb: %zu, end comb: %zu\n", rankn, start_combination, end_combination);
#endif
          int64_t counter = end_combination;
          for(int j=4; j>0; j--) {
            uint64_t n_comb = get_num_combinations(n_neighbors, j);
            end_comb_subgraph[j] = n_comb;
            counter -= n_comb;
            if(counter > start_combination) {
              start_comb_subgraph[j] = 0;
            } else {
              start_comb_subgraph[j] = start_combination-counter;
              break;
            }
          }
#ifdef DEBUG
printf("Rank %d: [start,end): ", rankn);
for(int j=0; j<5; j++) {
  printf("[%zu,%zu) ", start_comb_subgraph[j], end_comb_subgraph[j]);
}
printf("\n");
#endif
        } else if(first_chunk) {
          // First chunk
          start_combination = 0;
          end_combination = k_interval;
#ifdef DEBUG
printf("Rank %d: First chunk: start comb: %zu, end comb: %zu\n", rankn, start_combination, end_combination);
#endif
          uint64_t counter = k_interval;
          for(int j=1; j<5; j++) {
            uint64_t n_comb = get_num_combinations(n_neighbors, j);
            if(counter > n_comb) {
              end_comb_subgraph[j] = n_comb;
              counter -= n_comb;
            } else {
              end_comb_subgraph[j] = counter;
              break;
            }
          }
#ifdef DEBUG
printf("Rank %d: [start,end): ", rankn);
for(int j=0; j<5; j++) {
  printf("[%zu,%zu) ", start_comb_subgraph[j], end_comb_subgraph[j]);
}
printf("\n");
#endif
        } else if(middle_chunk) {
          // Middle chunk
          start_combination = start_chunk*k_interval;
          end_combination = start_combination+k_interval;
#ifdef DEBUG
printf("Rank %d: Middle chunk: start comb: %zu, end comb: %zu\n", rankn, start_combination, end_combination);
#endif
          uint64_t counter = 0;
          for(int j=1; j<5; j++) {
            uint64_t n_comb = get_num_combinations(n_neighbors, j);
#ifdef DEBUG
printf("Rank %d: Node count: %d, Number of combinations: %lu\n", rankn, j+1, n_comb);
#endif
            if(start_combination > counter && start_combination < counter+n_comb) {
              start_comb_subgraph[j] = start_combination-counter;
            }
            if(end_combination >= counter+n_comb && counter+n_comb >= start_combination) {
              end_comb_subgraph[j] = n_comb;
            }
            if(end_combination > counter && end_combination < counter+n_comb) {
              end_comb_subgraph[j] = end_combination-counter;
              break;
            }
            counter += n_comb;
          }
#ifdef DEBUG
printf("Rank %d: [start,end): ", rankn);
for(int j=0; j<5; j++) {
  printf("[%zu,%zu) ", start_comb_subgraph[j], end_comb_subgraph[j]);
}
printf("\n");
#endif
        }
        
        int node = curr_node;
        uint64_t num_bytes = 7*sizeof(scratch_int_view) + sizeof(scratch_subgraph_view) + 
                             sizeof(scratch_bool_view) + sizeof(bool)*5 + 
                             (max_neighbors+61)*sizeof(int);
#ifdef DEBUG
printf("Rank: %d: Per thread scratch %ld\n", rankn, num_bytes);
#endif
        Kokkos::TeamPolicy<> team_policy = Kokkos::TeamPolicy<>(1, NUM_THREADS).set_scratch_size(1, Kokkos::PerThread((num_bytes)));
        for(int node_count = 1; node_count < 5; node_count++) {
          uint64_t s_comb = start_comb_subgraph[node_count];
          uint64_t e_comb = end_comb_subgraph[node_count];
#ifdef DEBUG
printf("Rank %d: Node count: %d, range: [%zu,%zu)\n", rankn, node_count+1, s_comb, e_comb);
#endif
          if(e_comb-s_comb > 0) {
            Kokkos::parallel_for("Calculate GDV", team_policy, KOKKOS_LAMBDA(member_type team_member) {
              scratch_int_view scratch_view(team_member.thread_scratch(1), max_neighbors);
              scratch_int_view indices_subview(team_member.thread_scratch(1), 5);
              scratch_int_view combination_subview(team_member.thread_scratch(1), 5);
              scratch_int_view sgraph_dist_subview(team_member.thread_scratch(1), 6);
              scratch_int_view sgraph_deg_subview(team_member.thread_scratch(1), 5);
              scratch_int_view queue_subview(team_member.thread_scratch(1), 5);
              scratch_int_view distance_subview(team_member.thread_scratch(1), 5);
              scratch_subgraph_view subgraph_subview(team_member.thread_scratch(1), 5, 5);
              scratch_bool_view visited_subview(team_member.thread_scratch(1), 5);
              Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, e_comb-s_comb), [=] (uint64_t idx) {
                uint64_t combination_num = idx+s_comb;
//printf("Rank %d: node count: %d, Created combination number: %lu\n", rankn, node_count, combination_num);
                combination_from_position(scratch_view, combination_num, n_neighbors, node_count);
                for(int j=0; j<node_count; j++) {
                  combination_subview(j) = neighbor_scratch(scratch_view(j));
                }
                combination_subview(node_count) = node;
//printf("Rank %d: node count: %d, Created combination number: %lu\n", rankn, node_count, combination_num);
                kokkos_calculate_GDV2(team_member, node, node_count, graph, orbits, num_neighbors(curr_node), combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, graph_GDV);
//printf("Rank %d: node count: %d, update gdvs for combination number: %lu\n", rankn, node_count, combination_num);
              });
            });
          }
          Kokkos::fence();
//          printf("Rank %d: Done with subgraphs of size %d\n", rankn, node_count+1);
        }
#ifdef DEBUG
Kokkos::fence();
printf("Rank %d done with chunk %d\n", rankn, chunk_idx);
#endif
      }
#ifdef AUTO_CHECKPOINT
}, filt);
#endif
#ifdef DETAILED_TIMERS
      Kokkos::fence();
      chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
      chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2-t1);
//      printf("Rank: %d: Done with chunk: %u, time: %f\n", rankn, chunk_idx, time_span.count());
#endif
      double curr_time = MPI_Wtime();

      if((offset == intervals_per_rank-1) || (curr_time - prior_time > TIME_INTERVAL) || (TIME_INTERVAL == 0)) {
        std::string filename = label + "." + std::to_string(chkpt_counter) + ".hashtree.incr_chkpt";
        std::string logname = label + "." + std::to_string(chkpt_counter);
        uint64_t gdv_len = graph_GDV.span()*sizeof(uint32_t);
        Kokkos::View<uint8_t*>::HostMirror diff_h;
        auto gdv_h = Kokkos::create_mirror_view(graph_GDV);
        Kokkos::deep_copy(gdv_h, graph_GDV);
        std::string ref_digest = calculate_digest_host(gdv_h);

        deduplicator.checkpoint(dedup_mode, (uint8_t*)(graph_GDV.data()), gdv_len, diff_h, logname, chkpt_counter==0);
        Kokkos::fence();
        diffs.push_back(diff_h);
        deduplicator.restart(dedup_mode, (uint8_t*)(graph_GDV.data()), gdv_len, diffs, logname, chkpt_counter);

        Kokkos::deep_copy(gdv_h, graph_GDV);
        std::string digest = calculate_digest_host(gdv_h);

        if(ref_digest != digest)
          std::cout << "Rank " << rankn << ", Chunk " << chkpt_counter << ": Restart failed! Mismatch digests " << ref_digest << " vs " << digest << std::endl;

        printf("Rank: %d: Done with chunk: %u, chkpt counter %u, time: %f\n", rankn, offset, chkpt_counter, curr_time-prior_time);
        prior_time = curr_time;
        chkpt_counter += 1;
      }
      Kokkos::fence();

      offset++;
    }
    auto graph_GDV_h = Kokkos::create_mirror_view(graph_GDV);
    Kokkos::deep_copy(graph_GDV_h, graph_GDV);
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG
printf("Trying to reduce GDVs\n");
#endif
    if(rankn == 0) {
      MPI_Reduce(MPI_IN_PLACE, graph_GDV_h.data(), graph_GDV.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
      MPI_Reduce(graph_GDV_h.data(), graph_GDV_h.data(), graph_GDV.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    Kokkos::deep_copy(graph_GDV, graph_GDV_h);
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG
printf("Reduced GDVs\n");
#endif
  MPI_Barrier(MPI_COMM_WORLD);

  vec_calc_communication_time[graph_counter - 1] = process_ends_communication - vec_calc_start;

  #ifdef DEBUG
    cout << "Finished GDV Vector Calc on Rank: " << rankn << endl;
  #endif
}


//void 
//kokkos_GDV_vector_calculation_messy(const matrix_type& graph, 
//                              GDVs& graph_GDV, 
//                              const Orbits& orbits, 
//                              const char* graph_name, 
//                              int graph_counter) {
//  // Set up MPI parallelization               
//  int comm_size, rankn;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rankn);
//  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
////  int tag = 11;
//  uint32_t graph_size = graph.numRows();
//  std::string label(graph_name);
//  label += "_Rank";
//  label += std::to_string(rankn);
//
//#ifdef DIRTY_PAGE_TRACKING
//  pid_t pid = getpid();
//  uint64_t page_size = sysconf(_SC_PAGE_SIZE);
//#endif
//
//#ifdef DEBUG
//  cout << "Starting GDV vector calculation\n";
//  cout << "# of processes : " << comm_size << endl;
//#endif
//
//  double process_ends_communication;
//  double vec_calc_start = MPI_Wtime();
//
//  Kokkos::View<uint64_t*> num_combinations("Number of combinations", graph.numRows());
//  Kokkos::View<uint32_t*> num_neighbors("Number of neighbors", graph.numRows());
//  uint64_t k_interval = CHECKPOINT_INTERVAL;
//
//  typedef Kokkos::DefaultExecutionSpace::scratch_memory_space ScratchSpace;
//  typedef Kokkos::View<int*, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_int_view;
//  typedef Kokkos::View<bool*, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_bool_view;
////  int scratch_size = (4+sizeof(int)*4)*graph.numRows();
////printf("Scratch size: %d\n", scratch_size);
////  Kokkos::TeamPolicy<> team_policy = Kokkos::TeamPolicy<>(1, NUM_THREADS).set_scratch_size(1, Kokkos::PerThread((4+sizeof(int)*sizeof(bool))*graph.numRows()));
//  Kokkos::TeamPolicy<> team_policy = Kokkos::TeamPolicy<>(1, NUM_THREADS);
//  using team_member_type = Kokkos::TeamPolicy<>::member_type;
//#ifdef DEBUG
//  printf("Number of threads: %d\n", team_policy.team_size());
//#endif
////  using team_member_type = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Dynamic>>::member_type;
////  printf("Num threads: %d, num rows: %d\n", team_policy.team_size(), graph.numRows());
////
//  Kokkos::View<int**> neighbor_scratch("Neighbors", team_policy.team_size(), graph.numRows());
//  Kokkos::View<bool**> visited_scratch("Visited", team_policy.team_size(), graph.numRows());
//  Kokkos::deep_copy(visited_scratch, false);
//#ifdef DEBUG
//  printf("Starting combination counter kernel\n");
//#endif
//
//  uint32_t vertices_per_proc = graph.numRows()/comm_size;
//  uint64_t start_offset = rankn*vertices_per_proc;
//  uint64_t end_offset = start_offset+vertices_per_proc;
//  if(end_offset > graph.numRows())
//    end_offset = graph.numRows();
//  int* recv_count = new int[comm_size];
//  int* displs = new int[comm_size];
//  int running_count = 0;
//  for(uint64_t i=0; i<comm_size-1; i++) {
//    recv_count[i] = vertices_per_proc;
//    displs[i] = running_count;
//    running_count += vertices_per_proc;
//  }
//  recv_count[uint64_t (comm_size-1)] = graph.numRows()-(comm_size-1)*vertices_per_proc;
//  displs[comm_size-1] = running_count;
//#ifdef DEBUG
//  for(int i=0; i<comm_size; i++) {
//    printf("(%d,%d) ", recv_count[i], displs[i]);
//  }
//  printf("\n");
//  printf("(%u,%u)\n", start_offset, end_offset);
//#endif
//
//  chrono::steady_clock::time_point count_comb_beg = chrono::steady_clock::now();
//  Kokkos::parallel_for("Get # of combinations", team_policy, KOKKOS_LAMBDA(team_member_type team_member) {
////    scratch_int_view neighbor_subview(team_member.thread_scratch(1), graph.numRows());
////    scratch_bool_view visited_subview(team_member.thread_scratch(1), graph.numRows());
////    for(int j=0; j<graph.numRows(); j++) {
////      visited_subview(j) = false;
////    }
//    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, end_offset-start_offset), [=] (int n) {
//      uint64_t node = start_offset+n;
//      auto neighbor_subview = Kokkos::subview(neighbor_scratch, team_member.team_rank(), Kokkos::ALL());
//      auto visited_subview = Kokkos::subview(visited_scratch, team_member.team_rank(), Kokkos::ALL());
//      uint32_t n_neighbors = 0;
////      n_neighbors = EssensKokkos::get_num_neighbors(graph, node, 4, neighbor_subview);
//      n_neighbors = EssensKokkos::get_num_neighbors(graph, node, 4, neighbor_subview, visited_subview);
////      n_neighbors = EssensKokkos::get_num_neighbors_visited(graph, node, 4, visited_subview);
//      num_neighbors(node) = n_neighbors;
//      for(int i=1; i<5; i++) {
//        num_combinations(node) += get_num_combinations(n_neighbors, i);
//      }
//    });
//  });
//  Kokkos::fence();
//  chrono::steady_clock::time_point count_comb_end = chrono::steady_clock::now();
//  chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(count_comb_end-count_comb_beg);
//  printf("Time spent counting combinations: %f\n", time_span.count());
//
//
//#ifdef DEBUG
//  printf("Counted number of combinations for each vertex\n");
//#endif
//
//  // Find the maximum number of neighbors within distance 4
//  uint32_t max_neighbors = 0;
//  Kokkos::parallel_reduce("Find max neighbors", graph.numRows(), KOKKOS_LAMBDA(const uint32_t& x, uint32_t& deg) {
//    if(num_neighbors(x) > deg)
//      deg = num_neighbors(x);
//  }, Kokkos::Max<uint32_t>(max_neighbors));
//  printf("Max degree: %d\n", max_neighbors);
//
//  // Get the total # of combinations of potential subgraphs
//  uint64_t total_combinations = 0;
//  uint64_t local_combinations = 0;
//  Kokkos::parallel_reduce("Number of combinations", end_offset-start_offset, 
//  KOKKOS_LAMBDA(const int i, uint64_t& update) {
//    update += num_combinations(i+start_offset);
//  }, local_combinations);
//  Kokkos::View<uint64_t*> send_comb("Send buffer", end_offset-start_offset);
//  Kokkos::View<uint64_t*>::HostMirror send_comb_host = Kokkos::create_mirror_view(send_comb);
//  auto s_view = Kokkos::subview(num_combinations, Kokkos::make_pair(start_offset, end_offset));
//  Kokkos::deep_copy(send_comb, s_view);
//  Kokkos::deep_copy(send_comb_host, send_comb);
//  Kokkos::View<uint64_t*>::HostMirror num_combinations_host = Kokkos::create_mirror_view(num_combinations);
//  Kokkos::deep_copy(num_combinations_host, num_combinations);
//  MPI_Barrier(MPI_COMM_WORLD);
//  MPI_Allgatherv(send_comb_host.data(), 
//                  end_offset-start_offset, 
//                  MPI_UNSIGNED_LONG, 
//                  num_combinations_host.data(), 
//                  recv_count, 
//                  displs, 
//                  MPI_UNSIGNED_LONG, 
//                  MPI_COMM_WORLD);
//  delete[] recv_count, displs;
//  MPI_Allreduce(&local_combinations, &total_combinations, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
//  MPI_Barrier(MPI_COMM_WORLD);
//  Kokkos::deep_copy(num_combinations, num_combinations_host);
//
//  // Split combinations into intervals
//  uint64_t num_intervals = total_combinations/k_interval;
//  if(num_intervals*k_interval < total_combinations) {
//    num_intervals++;
//  }
//  if(rankn == 0) {
//    printf("Computed # of combinations: (%lu) split into %lu groups\n", total_combinations, num_intervals);
//  }
////  Kokkos::View<uint32_t*> starts("Start indices", num_intervals+1);
////  Kokkos::View<uint32_t*>::HostMirror starts_host = Kokkos::create_mirror_view(starts);
////  Kokkos::deep_copy(starts, 0);
////  Kokkos::deep_copy(starts_host, 0);
////  Kokkos::parallel_for("Find starts", Kokkos::RangePolicy<>(0, num_intervals+1), KOKKOS_LAMBDA(const uint64_t i) {
////    if(i == num_intervals) {
////      starts(i) = graph.numRows();
////    } else {
////      uint64_t threshold = i*k_interval;
////      uint64_t counter = 0;
////      for(uint32_t j=0; j<num_combinations.extent(0); j++) {
////        counter += num_combinations(j);
////        if(counter > threshold) {
////          starts(i) = j;
////          break;
////        }
////      }
////    }
////  });
////  Kokkos::fence();
////  Kokkos::deep_copy(starts_host, starts);
//  Kokkos::fence();
//#ifdef DEBUG
//  printf("Rank %d: ", rankn);
//  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const uint32_t _i) {
//    for(uint32_t i=0; i<num_intervals+1; i++) {
////      printf("%ld ", starts_host(i));
////      printf("%ld ", get_node_from_combinations(graph, num_combinations_host, k_interval, i));
//    }
//  });
//  printf("\n");
//#endif
//
//  // Calculate GDVs
//  if(comm_size == 1)
//  {
//    using member_type = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Dynamic>>::member_type;
//
////    int i=0;
//#ifdef AUTO_CHECKPOINT
//    auto ctx = KokkosResilience::make_context(MPI_COMM_SELF, "/home/ntan1/Src_Fido_Kokkos/fido.json");
//    printf("Created context\n");
//    const auto filt = KokkosResilience::Filter::NthIterationFilter(1);
//    printf("Created filter\n");
//    i = KokkosResilience::latest_version(*ctx, graph_name);
//    if(i < 0)
//      i = 0;
//    printf("Got latest counter %d\n", i);
//#endif
//#ifdef ANALYSIS
//    int checkpoint_num = 0;
//#ifdef HASH_ANALYSIS
////    uint64_t page_size = sysconf(_SC_PAGE_SIZE);
//    uint64_t page_size = CHUNK_SIZE;
//    Kokkos::View<uint64_t*> old_hashes("Old hashes", graph_GDV.span()*sizeof(uint32_t)/page_size);
//    Kokkos::View<uint64_t*> new_hashes("New hashes", graph_GDV.span()*sizeof(uint32_t)/page_size);
//    auto& old_data = old_hashes;
//    auto& new_data = new_hashes;
//#endif
//#ifdef SCAN_ANALYSIS
//    GDVs old_updates("Old updates", graph_GDV.layout());
//    auto& old_data = old_updates;
//    auto& new_data = graph_GDV;
//#endif
//#endif
//
//#ifdef DEDUPLICATION
//    std::fstream fs;
//    HASH_FUNC hash_func;
//    using hash_digest_type = HASH_FUNC::DIGEST_TYPE;
//    Kokkos::View<hash_digest_type**, Kokkos::DefaultHostExecutionSpace> hist_hashes;
//    Kokkos::View<hash_digest_type*> prev_hashes, curr_hashes;
//    Kokkos::View<hash_digest_type*, Kokkos::DefaultHostExecutionSpace> prev_hashes_h, curr_hashes_h;
//    Kokkos::View<uint64_t*> changed_regions;
//    Kokkos::View<uint64_t*, Kokkos::DefaultHostExecutionSpace> changed_regions_h;
//    int chunk_size = CHUNK_SIZE;
//    printf("Chunk size: %d\n", chunk_size);
//    uint64_t num_hashes = graph_GDV.span()*sizeof(uint32_t) / chunk_size;
//    if(chunk_size*num_hashes < graph_GDV.span()*sizeof(uint32_t))
//      num_hashes += 1;
//    uint64_t num_elements = num_hashes*hash_func.digest_size()/sizeof(hash_digest_type);
//    printf("Number of hashes: %zd\n", num_hashes);
//    printf("Number of elements in hash view: %zd\n", num_elements);
//    // Storage for hashes
//printf("Hash buff size: %u\n", num_hashes*20);
//#ifdef DEDUP_CASE0
//    changed_regions_h = Kokkos::View<uint64_t*, Kokkos::DefaultHostExecutionSpace>("Changed regions host", num_hashes);
//    std::string logfile = std::string("case0.chunk_size_") + 
//                          std::to_string(CHUNK_SIZE) + 
//                          std::string(".") + 
//                          std::string(graph_name);
//    fs.open(logfile+std::string(".csv"), std::fstream::out | std::fstream::app);
//    printf("Deduplication on CPU (case 0)\n");
//    fs << "Chunk, Number of changed regions, Copy GDVs to CPU, Find differences, Create diff, Total time\n";
//    fs.close();
//#elif defined(DEDUP_CASE1)
//    hist_hashes = Kokkos::View<hash_digest_type**, Kokkos::DefaultHostExecutionSpace>("History of hashes", NUM_CHECKPOINTS, num_elements);
//    prev_hashes_h = Kokkos::View<hash_digest_type*, Kokkos::DefaultHostExecutionSpace>("Previous hashes host", num_elements);
//    curr_hashes_h = Kokkos::View<hash_digest_type*, Kokkos::DefaultHostExecutionSpace>("Current hashes host", num_elements);
//    changed_regions_h = Kokkos::View<uint64_t*, Kokkos::DefaultHostExecutionSpace>("Changed regions host", num_hashes);
//    std::string logfile = std::string("case1.chunk_size_") + 
//                          std::to_string(CHUNK_SIZE) + 
//                          std::string(".") + 
//                          hash_func.hash_name() + 
//                          std::string(".") + 
//                          std::string(graph_name);
//    fs.open(logfile+std::string(".csv"), std::fstream::out | std::fstream::app);
//    printf("Deduplication on CPU (case 1)\n");
//    fs << "Chunk, Number of changed regions, Copy GDVs to CPU, Compute hashes, Compare hashes, Create diff, Total time\n";
//    fs.close();
//#elif defined(DEDUP_CASE1_ASYNC)
//    prev_hashes_h = Kokkos::View<hash_digest_type*, Kokkos::DefaultHostExecutionSpace>("Previous hashes host", num_elements);
//    curr_hashes_h = Kokkos::View<hash_digest_type*, Kokkos::DefaultHostExecutionSpace>("Current hashes host", num_elements);
//    changed_regions_h = Kokkos::View<uint64_t*, Kokkos::DefaultHostExecutionSpace>("Changed regions host", num_hashes);
//    std::string logfile = std::string("case1_async.chunk_size_") + 
//                          std::to_string(CHUNK_SIZE) + 
//                          std::string(".") + 
//                          hash_func.hash_name() + 
//                          std::string(".") + 
//                          std::string(graph_name);
//    fs.open(logfile+std::string(".csv"), std::fstream::out | std::fstream::app);
//    printf("Deduplication on CPU (case 1)\n");
//    fs << "Chunk, Number of changed regions, Copy GDVs to CPU, Compute hashes, Compare hashes, Create diff, Total time\n";
//    fs.close();
//#elif defined(DEDUP_CASE2)
//    hist_hashes = Kokkos::View<hash_digest_type**, Kokkos::DefaultHostExecutionSpace>("History of hashes", NUM_CHECKPOINTS, num_elements);
//    prev_hashes = Kokkos::View<hash_digest_type*>("Previous hashes", num_elements);
//    curr_hashes = Kokkos::View<hash_digest_type*>("Current hashes", num_elements);
//    Kokkos::deep_copy(prev_hashes, 0);
//    Kokkos::deep_copy(curr_hashes, 0);
//    prev_hashes_h = Kokkos::View<hash_digest_type*, Kokkos::DefaultHostExecutionSpace>("Previous hashes host", num_elements);
//    changed_regions = Kokkos::View<uint64_t*>("Regions that changed", num_hashes);
//    changed_regions_h = Kokkos::View<uint64_t*, Kokkos::DefaultHostExecutionSpace>("Changed regions host", num_hashes);
//    std::string logfile = std::string("case2.chunk_size_") + 
//                          std::to_string(CHUNK_SIZE) + 
//                          std::string(".") + 
//                          hash_func.hash_name() + 
//                          std::string(".") + 
//                          std::string(graph_name);
//    fs.open(logfile+std::string(".csv"), std::fstream::out | std::fstream::app);
//    printf("Deduplication on GPU (case 2)\n");
//    fs << "Chunk, Number of changed regions, Copy hashes to GPU, Compute hashes, Compare hashes, Copy changes to CPU, Create diff, Copy diff to CPU, Total time\n";
//    fs.close();
//#elif defined(DEDUP_CASE3)
//    curr_hashes = Kokkos::View<hash_digest_type*>("Current hashes", num_elements);
//    Kokkos::deep_copy(prev_hashes, 0);
//    Kokkos::deep_copy(curr_hashes, 0);
//    prev_hashes_h = Kokkos::View<hash_digest_type*, Kokkos::DefaultHostExecutionSpace>("Previous hashes host", num_elements);
//    curr_hashes_h = Kokkos::View<hash_digest_type*, Kokkos::DefaultHostExecutionSpace>("Current hashes host", num_elements);
//    changed_regions = Kokkos::View<uint64_t*>("Regions that changed", num_hashes);
//    changed_regions_h = Kokkos::View<uint64_t*, Kokkos::DefaultHostExecutionSpace>("Changed regions host", num_hashes);
//    std::string logfile = std::string("case3.chunk_size_") + 
//                          std::to_string(CHUNK_SIZE) + 
//                          std::string(".") + 
//                          hash_func.hash_name() + 
//                          std::string(".") + 
//                          std::string(graph_name);
//    fs.open(logfile+std::string(".csv"), std::fstream::out | std::fstream::app);
//    printf("Deduplication on GPU (case 3)\n");
//    fs << "Chunk, Number of changed regions, Compute hashes, Copy hashes to CPU, Compare hashes, Create diff, Copy diff to CPU, Total time\n";
//    fs.close();
//#endif
//    using namespace std::chrono;
////    prev_hashes = Kokkos::View<uint8_t*>("Previous hashes", num_hashes*20);
////    curr_hashes = Kokkos::View<uint8_t*>("Current hashes", num_hashes*20);
////    Kokkos::deep_copy(prev_hashes, 0);
////    Kokkos::deep_copy(curr_hashes, 0);
////    prev_hashes_h = Kokkos::View<uint8_t*, Kokkos::DefaultHostExecutionSpace>("Previous hashes host", num_hashes*20);
////    curr_hashes_h = Kokkos::View<uint8_t*, Kokkos::DefaultHostExecutionSpace>("Current hashes host", num_hashes*20);
////    changed_regions = Kokkos::View<uint64_t*>("Regions that changed", num_hashes);
////    changed_regions_h = Kokkos::View<uint64_t*, Kokkos::DefaultHostExecutionSpace>(changed_regions);
//    GDVs::HostMirror current_gdvs_h = Kokkos::create_mirror_view(graph_GDV);
//    Kokkos::View<uint8_t*, Kokkos::HostSpace> diff_buff_h("Diff", graph_GDV.span()*sizeof(uint32_t));
//    Kokkos::fence();
////printf("Allocated memory for deduplication\n");
//#endif
//
//int num_leagues = 1;
//      Kokkos::TeamPolicy<> bundle_policy = Kokkos::TeamPolicy<>(num_leagues, NUM_THREADS);
//      Kokkos::View<int**> neighbor_view("Neighbor scratch", bundle_policy.team_size(), max_neighbors);
//      Kokkos::View<int**> indices_view("Indices scratch", bundle_policy.team_size(), 5);
//      Kokkos::View<int**> combination_view("Combination scratch", bundle_policy.team_size(), 5);
//      Kokkos::View<int**> sgraph_dist_view("Subgraph Distance scratch", bundle_policy.team_size(), orbits.distance.extent(0));
//      Kokkos::View<int**> sgraph_deg_view("Degree scratch", bundle_policy.team_size(), 5);
//      Kokkos::View<int**> queue_view("Queue scratch", bundle_policy.team_size(), 5);
//      Kokkos::View<int**> distance_view("Distance scratch", bundle_policy.team_size(), 5);
//      Kokkos::View<int***> subgraph_view("Subgraph scratch", bundle_policy.team_size(), 5, 5);
//      Kokkos::View<bool**> visited_view("Visited scratch", bundle_policy.team_size(), 5);
//
//    for(int i=0; i<num_intervals; i++) {
//
//#ifdef DIRTY_PAGE_TRACKING
//      printf("Chunk %d\n", i);
//      reset_dirty_bit(pid);
//#endif
//
//#ifdef AUTO_CHECKPOINT
//      KokkosResilience::checkpoint(*ctx, graph_name, i, [=] () mutable {
//#endif
//#ifdef DETAILED_TIMERS
//      chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
//#endif
//
//      // Divide work into different leagues
////      const int num_leagues = 1;
//      int start_node = get_node_from_combinations(graph, num_combinations_host, k_interval, i);
//      int end_node = 0;
//      if((i+1)*k_interval > total_combinations) {
//        end_node = graph.numRows();
//      } else {
//        end_node = get_node_from_combinations(graph, num_combinations_host, k_interval, (i+1));
//      }
//      int per_league = (end_node - start_node)/num_leagues;
//      if(per_league*num_leagues < end_node-start_node)
//        per_league += 1;
////      int per_league = (starts_host(i+1)-starts_host(i))/num_leagues;
////      if(per_league*num_leagues < starts_host(i+1)-starts_host(i))
////        per_league += 1;
////#ifdef DEBUG
////      printf("start: %d, end: %d\n", starts_host(i), starts_host(i+1));
//      printf("start: %d, end: %d\n", start_node, end_node);
//      printf("Num leagues: %d, per league: %d\n", num_leagues, per_league);
////#endif
//
//      typedef Kokkos::DefaultExecutionSpace::scratch_memory_space ScratchSpace;
//      typedef Kokkos::View<int*, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_int_view;
//      typedef Kokkos::View<int**, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_subgraph_view;
//      typedef Kokkos::View<bool*, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_bool_view;
//
////      int num_bytes = (8+max_neighbors+55+orbits.distance.span())*sizeof(int);
////      Kokkos::TeamPolicy<> bundle_policy = Kokkos::TeamPolicy<>(num_leagues, NUM_THREADS).set_scratch_size(1, Kokkos::PerThread(((8+max_neighbors+55+orbits.distance.span())*sizeof(int))));
//      Kokkos::TeamPolicy<> bundle_policy = Kokkos::TeamPolicy<>(num_leagues, NUM_THREADS);
//
//      Kokkos::parallel_for("Calculate GDV bundle", bundle_policy, KOKKOS_LAMBDA(member_type team_member) {
//
//        auto neighbor_subview    = Kokkos::subview(neighbor_view,    team_member.team_rank(), Kokkos::ALL());
//        auto indices_subview     = Kokkos::subview(indices_view,     team_member.team_rank(), Kokkos::ALL());
//        auto combination_subview = Kokkos::subview(combination_view, team_member.team_rank(), Kokkos::ALL());
//        auto sgraph_dist_subview = Kokkos::subview(sgraph_dist_view, team_member.team_rank(), Kokkos::ALL());
//        auto sgraph_deg_subview  = Kokkos::subview(sgraph_deg_view,  team_member.team_rank(), Kokkos::ALL());
//        auto queue_subview       = Kokkos::subview(queue_view,       team_member.team_rank(), Kokkos::ALL());
//        auto distance_subview    = Kokkos::subview(distance_view,    team_member.team_rank(), Kokkos::ALL());
//        auto subgraph_subview    = Kokkos::subview(subgraph_view,    team_member.team_rank(), Kokkos::ALL(), Kokkos::ALL());
//        auto visited_subview     = Kokkos::subview(visited_view,     team_member.team_rank(), Kokkos::ALL());
//
////        scratch_int_view neighbor_subview(team_member.thread_scratch(1), max_neighbors);
////        scratch_int_view indices_subview(team_member.thread_scratch(1), 5);
////        scratch_int_view combination_subview(team_member.thread_scratch(1), 5);
////        scratch_int_view sgraph_dist_subview(team_member.thread_scratch(1), orbits.distance.extent(0));
////        scratch_int_view sgraph_deg_subview(team_member.thread_scratch(1), 5);
////        scratch_int_view queue_subview(team_member.thread_scratch(1), 5);
////        scratch_int_view distance_subview(team_member.thread_scratch(1), 5);
////        scratch_subgraph_view subgraph_subview(team_member.thread_scratch(1), 5, 5);
////        scratch_bool_view visited_subview(team_member.thread_scratch(1), 5);
//
//        Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, per_league), [=] (int n_) {
//          int n_offset = team_member.league_rank()*per_league + n_;
////          int node = starts(i) + n_offset;
////          if(node < starts(i+1)) {
//          int node = start_node + n_offset;
//          if(node < end_node) {
//            kokkos_calculate_GDV( team_member, 
//                                  node, 
//                                  graph, 
//                                  orbits, 
//                                  neighbor_subview, 
//                                  indices_subview, 
//                                  combination_subview, 
//                                  sgraph_dist_subview, 
//                                  sgraph_deg_subview, 
//                                  subgraph_subview, 
//                                  visited_subview, 
//                                  queue_subview, 
//                                  distance_subview, 
//                                  graph_GDV);
////                                  metrics_sa);
//          }
//        });
//      });
//#ifndef DEDUP_CASE1_ASYNC
//      Kokkos::fence();
//#endif
////      Kokkos::Experimental::contribute(graph_GDV, metrics_sa);
////      metrics_sa.reset_except(graph_GDV);
//#ifndef DEDUP_CASE1_ASYNC
//      Kokkos::fence();
//#endif
//
//#ifdef ANALYSIS
//std::string region_log("region-data-k-");
//region_log = region_log + std::to_string(CHECKPOINT_INTERVAL) + "-blocksize-" + std::to_string(CHUNK_SIZE) + std::string("-graph-");
//region_log = region_log + std::string(graph_name) + std::string(".log");
//std::fstream fs(region_log, std::fstream::out|std::fstream::app);
//// Hash based change tracking
//uint32_t num_changed = 0;
//#ifdef HASH_ANALYSIS
////auto test = make_layout_left(graph_GDV);
////generate_hashes(test, new_hashes, page_size);
//generate_hashes(graph_GDV, new_hashes, page_size);
//#endif
//num_changed = print_changed_blocks(fs, new_data, old_data);
//auto contiguous_regions = print_contiguous_regions(fs, new_data, old_data);
//
//#ifdef HASH_ANALYSIS
//Kokkos::deep_copy(old_data, new_data);
//Kokkos::deep_copy(new_data, 0);
//#endif
//#ifdef SCAN_ANALYSIS
//Kokkos::deep_copy(old_data, new_data);
//#endif
//#endif
//
//#ifdef DIRTY_PAGE_TRACKING
//      uintptr_t start_addr = reinterpret_cast<uintptr_t>(graph_GDV.data());
//      uintptr_t end_addr = reinterpret_cast<uintptr_t>(graph_GDV.data() + graph_GDV.span());
//      bool range_dirty = address_range_dirty(pid, start_addr, end_addr);
//      if(range_dirty) {
//        printf("Detected changes\n");
//      }
//      uint64_t num_pages = ((end_addr-start_addr)/page_size);
//      uint64_t* page_list = (uint64_t*) malloc(sizeof(uint64_t)*num_pages);
//      for(int j=0; j<num_pages; j++) {
//        page_list[j] = 0; 
//      }
//      bool dirty = get_dirty_pages(pid, start_addr, end_addr, page_list);
//      printf("List of dirty pages: ");
//      for(int j=0; j<num_pages; j++) {
//        if(page_list[j] > 0) {
//          printf("%llu, ", page_list[j]);
//        }
//      }
//      printf("\n");
//#endif
//
//#ifdef DETAILED_TIMERS
//      Kokkos::fence();
//      chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
//      chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2-t1);
//      printf("Done with chunk: %d, time: %f\n", i, time_span.count());
//#endif
//
//#ifdef AUTO_CHECKPOINT
//}, filt);
//#endif
//
//#ifdef DEDUPLICATION
//#ifdef DEDUP_CASE0
//  dedup_case0(graph_GDV, current_gdvs_h,
//              changed_regions, changed_regions_h,
//              diff_buff_h, logfile, chunk_size, num_hashes, i);
//#endif
//#ifdef DEDUP_CASE1
//  dedup_case1( 
//              graph_GDV, current_gdvs_h, 
//              prev_hashes, prev_hashes_h, 
//              curr_hashes, curr_hashes_h, 
//              changed_regions, changed_regions_h, 
//              diff_buff_h, logfile, chunk_size, num_hashes, i,
//              hist_hashes);
//#endif
//#ifdef DEDUP_CASE1_ASYNC
//  dedup_case1_async(graph_GDV, current_gdvs_h, 
//              prev_hashes, prev_hashes_h, 
//              curr_hashes, curr_hashes_h, 
//              changed_regions, changed_regions_h, 
//              diff_buff_h, logfile, chunk_size, num_hashes, num_intervals, i);
//#endif
//#ifdef DEDUP_CASE2
//  dedup_case2(graph_GDV, current_gdvs_h, 
//              prev_hashes, prev_hashes_h, 
//              curr_hashes, curr_hashes_h, 
//              changed_regions, changed_regions_h, 
//              diff_buff_h, logfile, chunk_size, num_hashes, i,
//              hist_hashes);
//#endif
//#ifdef DEDUP_CASE3
//  dedup_case3(graph_GDV, current_gdvs_h, 
//              prev_hashes, prev_hashes_h, 
//              curr_hashes, curr_hashes_h, 
//              changed_regions, changed_regions_h, 
//              diff_buff_h, logfile, chunk_size, num_hashes, i);
//#endif
//Kokkos::fence();
//#endif
//    }
//#ifdef DEDUPLICATION
//fs.close();
//#endif
//  }
//  else
//  {
//#ifdef ANALYSIS
//    int checkpoint_num = 0;
//#ifdef HASH_ANALYSIS
////    uint64_t page_size = sysconf(_SC_PAGE_SIZE);
//    uint64_t page_size = 4096;
//    Kokkos::View<uint64_t*> old_hashes("Old hashes", graph_GDV.span()*sizeof(uint32_t)/page_size);
//    Kokkos::View<uint64_t*> new_hashes("New hashes", graph_GDV.span()*sizeof(uint32_t)/page_size);
//    auto& old_data = old_hashes;
//    auto& new_data = new_hashes;
//#endif
//#ifdef SCAN_ANALYSIS
//    GDVs old_updates("Old updates", graph_GDV.layout());
//    auto& old_data = old_updates;
//    auto& new_data = graph_GDV;
//#endif
//#endif
//    Kokkos::deep_copy(graph_GDV, 0);
//    uint64_t intervals_per_rank = num_intervals/comm_size;
//    uint64_t remaining_intervals = num_intervals%comm_size;
//    if(rankn < remaining_intervals)
//      intervals_per_rank += 1;
//    uint64_t start_index = 0;
//    if(rankn < remaining_intervals) {
//      start_index = intervals_per_rank*rankn;
//    }
//    if(rankn >= remaining_intervals) {
//      start_index = ((intervals_per_rank+1)*remaining_intervals) + (intervals_per_rank*(rankn-remaining_intervals));
////      start_index = (intervals_per_rank+1)*remaining_intervals+intervals_per_rank*(rankn-remaining_intervals);
//    }
//    if(start_index+intervals_per_rank > num_intervals) {
//      intervals_per_rank = num_intervals - start_index;
//    }
//    int offset = 0;
//printf("Rank %d: start_index: %llu, intervals_per_rank: %llu\n", rankn, start_index, intervals_per_rank);
//    Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(1, NUM_THREADS);
//    using member_type = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type;
//
//    typedef Kokkos::DefaultExecutionSpace::scratch_memory_space ScratchSpace;
//    typedef Kokkos::View<int*, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_int_view;
//    typedef Kokkos::View<int**, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_subgraph_view;
//    typedef Kokkos::View<bool*, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> scratch_bool_view;
//
//    const int num_leagues = 1;
//    int num_bytes = 5*sizeof(bool) + (max_neighbors+max_neighbors+55+orbits.distance.span())*sizeof(int);
//    Kokkos::TeamPolicy<> team_policy = Kokkos::TeamPolicy<>(num_leagues, NUM_THREADS).set_scratch_size(0, Kokkos::PerThread((num_bytes)));
//
//#ifdef DEBUG
//printf("Rank %d allocated memory\n", rankn);
//#endif
//
//#ifdef AUTO_CHECKPOINT
//    auto ctx = KokkosResilience::make_context(MPI_COMM_SELF, "/home/ntan1/Src_Fido_Kokkos/fido.json");
//    printf("Created context\n");
//    const auto filt = KokkosResilience::Filter::NthIterationFilter(1);
//    printf("Created filter\n");
//    offset = KokkosResilience::latest_version(*ctx, label);
//    if(offset < 0)
//      offset = 0;
//    printf("Got latest counter %d\n", offset);
//#endif
//
//    while(offset < intervals_per_rank) {
//#ifdef DETAILED_TIMERS
//      chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
//#endif
//      uint32_t chunk_idx = start_index+offset;
//
//#ifdef AUTO_CHECKPOINT
//KokkosResilience::checkpoint(*ctx, label, offset, [=] () mutable {
//#endif
//#ifdef DIRTY_PAGE_TRACKING
//      reset_dirty_bit(pid);
//#endif
//      int prev_node  = get_node_from_combinations(graph, num_combinations_host, k_interval, (chunk_idx-1));
//      int curr_node = get_node_from_combinations(graph, num_combinations_host, k_interval, chunk_idx);
//      int next_node = 0;
//      if((chunk_idx+1)*k_interval > total_combinations) {
//        next_node = graph.numRows();
//      } else {
//        next_node = get_node_from_combinations(graph, num_combinations_host, k_interval, (chunk_idx+1));
//      }
//#ifdef DEBUG
//      printf("Rank %d: Prev node: %d, Curr node: %d, Next node: %d\n", rankn, prev_node, curr_node, next_node);
//#endif
//      bool multi_node = (chunk_idx+1 == num_intervals+1 || next_node != curr_node) && 
//                        (chunk_idx==0 || prev_node != curr_node);
////      bool multi_node = (chunk_idx+1 == num_intervals+1 || starts_host(chunk_idx+1) != starts_host(chunk_idx)) && 
////                        (chunk_idx==0 || starts_host(chunk_idx-1) != starts_host(chunk_idx));
//#ifdef DEBUG
//      printf("Rank %d: Multi node chunk?: %d\n", rankn, multi_node);
//#endif
//      if(multi_node) {
//printf("Multi node X\n");
////        int start_node = starts_host(chunk_idx);
//        int start_node = curr_node;
//        int end_node;
//        if(chunk_idx+1 == num_intervals+1) {
//          end_node = graph.numRows();
//        } else {
////          end_node = starts_host(chunk_idx+1);
//          end_node = next_node;
//        }
//
//        int per_league = (end_node-start_node)/num_leagues;
//        if(per_league*num_leagues < end_node-start_node)
//          per_league += 1;
//printf("Per league: %d\n", per_league);
//
//        Kokkos::parallel_for("Calculate GDV", team_policy, KOKKOS_LAMBDA(member_type team_member) {
//          scratch_int_view neighbor_subview(team_member.thread_scratch(0), max_neighbors);
//          scratch_int_view indices_subview(team_member.thread_scratch(0), 5);
//          scratch_int_view combination_subview(team_member.thread_scratch(0), 5);
//          scratch_int_view sgraph_dist_subview(team_member.thread_scratch(0), orbits.distance.extent(0));
//          scratch_int_view sgraph_deg_subview(team_member.thread_scratch(0), 5);
//          scratch_int_view queue_subview(team_member.thread_scratch(0), 5);
//          scratch_int_view distance_subview(team_member.thread_scratch(0), 5);
//          scratch_subgraph_view subgraph_subview(team_member.thread_scratch(0), 5, 5);
//          scratch_bool_view visited_subview(team_member.thread_scratch(0), 5);
//printf("Rank %d, start: %d, end: %d\n", rankn, team_member.league_rank()*per_league+0, team_member.league_rank()*per_league+per_league-1);
//          Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, per_league), [=] (int n) {
////            int node = start_node + n;
////            auto neighbor_subview = Kokkos::subview(all_neighbors, team_member.team_rank(), Kokkos::ALL());
////            auto indices_subview = Kokkos::subview(indices, team_member.team_rank(), Kokkos::ALL());
////            auto combination_subview = Kokkos::subview(combination_view, team_member.team_rank(), Kokkos::ALL());
////            auto sgraph_dist_subview = Kokkos::subview(sgraph_distance_signature, 
////                                                        team_member.team_rank(), 
////                                                        Kokkos::ALL());
////            auto sgraph_deg_subview = Kokkos::subview(sgraph_degree_signature, 
////                                                      team_member.team_rank(), 
////                                                      Kokkos::ALL());
////            auto subgraph_subview = Kokkos::subview(induced_subgraph, 
////                                                    team_member.team_rank(), 
////                                                    Kokkos::ALL(), 
////                                                    Kokkos::ALL());
////            auto visited_subview = Kokkos::subview(visited, team_member.team_rank(), Kokkos::ALL());
////            auto queue_subview = Kokkos::subview(queue, team_member.team_rank(), Kokkos::ALL());
////            auto distance_subview = Kokkos::subview(distance, team_member.team_rank(), Kokkos::ALL());
//            int n_offset = team_member.league_rank()*per_league + n;
//            int node = start_node + n_offset;
//printf("Rank %d, node: %d\n", rankn, node);
//            if(node < end_node) {
//              kokkos_calculate_GDV( team_member, 
//                                    node, 
//                                    graph, 
//                                    orbits, 
//                                    neighbor_subview, 
//                                    indices_subview, 
//                                    combination_subview, 
//                                    sgraph_dist_subview, 
//                                    sgraph_deg_subview, 
//                                    subgraph_subview, 
//                                    visited_subview, 
//                                    queue_subview, 
//                                    distance_subview, 
//                                    graph_GDV);
////                                    metrics_sa);
//            }
//          });
//        });
//#ifdef DEBUG
//printf("Rank %d done with chunk %d\n", rankn, chunk_idx);
//#endif
//      } else {
//printf("Split node\n");
//    Kokkos::View<int*> neighbor_scratch("Neighbor scratch", graph.numRows());
////        auto neighbor_scratch = Kokkos::subview(all_neighbors, 0, Kokkos::ALL());
////        int n_neighbors = EssensKokkos::find_neighbours(starts_host(chunk_idx), graph, 4, neighbor_scratch);
//        int n_neighbors = EssensKokkos::find_neighbours(curr_node, graph, 4, neighbor_scratch);
//
//        uint64_t start_combination;
//        uint64_t end_combination;
//        uint64_t start_comb_subgraph[5] = {0,0,0,0,0};
//        uint64_t end_comb_subgraph[5] = {0,0,0,0,0};
//        start_comb_subgraph[0] = 0;
//        end_comb_subgraph[0] = 0;
//        uint64_t start_chunk = 0;
//        for(int j=0; j<chunk_idx; j++) {
//          int u = get_node_from_combinations(graph, num_combinations_host, k_interval, j);
//          if(u == curr_node)
//            start_chunk++;
////          if(starts_host(j) == starts_host(chunk_idx))
//        }
//printf("Rank %d: start_chunk: %zu, k_interval: %zu\n", rankn, start_chunk, k_interval);
//        bool first_chunk=false, middle_chunk=false, last_chunk=false;
//        first_chunk = ((chunk_idx == 0) || 
//                            prev_node != curr_node) && 
//                            next_node == curr_node;
//        if(!first_chunk) {
//          last_chunk = prev_node == curr_node && 
//                       next_node != curr_node;
//          middle_chunk = prev_node == curr_node && 
//                         next_node == curr_node;
//        }
////        first_chunk = ((chunk_idx == 0) || 
////                            starts_host(chunk_idx-1) != starts_host(chunk_idx)) && 
////                            starts_host(chunk_idx+1) == starts_host(chunk_idx);
////        if(!first_chunk) {
////          last_chunk = starts_host(chunk_idx-1) == starts_host(chunk_idx) && 
////                       starts_host(chunk_idx+1) != starts_host(chunk_idx);
////          middle_chunk = starts_host(chunk_idx-1) == starts_host(chunk_idx) && 
////                         starts_host(chunk_idx+1) == starts_host(chunk_idx);
////        }
//        if(last_chunk) {
//          // Last chunk
//          start_combination = start_chunk*k_interval;
//          end_combination = num_combinations_host(curr_node);
//printf("Rank %d start combination: %zu, end combination: %zu\n", rankn, start_combination, end_combination);
////          end_combination = num_combinations_host(starts_host(chunk_idx));
//          int64_t counter = end_combination;
//          for(int j=4; j>0; j--) {
//            int64_t n_comb = get_num_combinations(n_neighbors, j);
//            end_comb_subgraph[j] = n_comb;
//            counter -= n_comb;
//            if(counter > start_combination) {
//              start_comb_subgraph[j] = 0;
//            } else {
//              start_comb_subgraph[j] = start_combination-counter;
//              break;
//            }
//          }
//printf("Rank %d: start combination: ", rankn);
//for(int j=0; j<5; j++) {
//  printf("%zu ", start_comb_subgraph[j]);
//}
//printf("\n");
//        } else if(first_chunk) {
//          // First chunk
//          start_combination = 0;
//          end_combination = k_interval;
//printf("Rank %d start combination: %zu, end combination: %zu\n", rankn, start_combination, end_combination);
//          uint64_t counter = k_interval;
//          for(int j=1; j<5; j++) {
//            int64_t n_comb = get_num_combinations(n_neighbors, j);
//            if(counter > n_comb) {
//              end_comb_subgraph[j] = n_comb;
//              counter -= n_comb;
//            } else {
//              end_comb_subgraph[j] = counter;
//              break;
//            }
//          }
//printf("Rank %d: start combination: ", rankn);
//for(int j=0; j<5; j++) {
//  printf("%zu ", start_comb_subgraph[j]);
//}
//printf("\n");
//        } else if(middle_chunk) {
//          // Middle chunk
//          start_combination = start_chunk*k_interval;
//          end_combination = start_combination+k_interval;
//printf("Rank %d start combination: %zu, end combination: %zu\n", rankn, start_combination, end_combination);
//          uint64_t counter = 0;
//          for(int j=1; j<5; j++) {
//            int64_t n_comb = get_num_combinations(n_neighbors, j);
//            if(start_combination > counter && start_combination < counter+n_comb) {
//              start_comb_subgraph[j] = start_combination-counter;
//            }
//            if(end_combination > counter+n_comb && counter+n_comb >= start_combination) {
//              end_comb_subgraph[j] = n_comb;
//            }
//            if(end_combination > counter && end_combination < counter+n_comb) {
//              end_comb_subgraph[j] = end_combination-counter;
//              break;
//            }
//            counter += n_comb;
//          }
//printf("Rank %d: start combination: ", rankn);
//for(int j=0; j<5; j++) {
//  printf("%zu ", start_comb_subgraph[j]);
//}
//printf("\n");
//        }
//printf("Prepared for kernel launch\n");        
////        int node = starts_host(chunk_idx);
//        int node = curr_node;
//        Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> team_policy(1, NUM_THREADS);
////        Kokkos::View<int**> scratch_block("Indices scratch", team_policy.team_size(), graph.numRows());
//        for(int node_count = 1; node_count < 5; node_count++) {
//          int64_t s_comb = start_comb_subgraph[node_count];
//          int64_t e_comb = end_comb_subgraph[node_count];
//          if(e_comb-s_comb > 0) {
//            Kokkos::parallel_for("Calculate GDV", team_policy, KOKKOS_LAMBDA(member_type team_member) {
//              Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, e_comb-s_comb), [=] (int64_t idx) {
//                int64_t combination_num = idx+s_comb;
////                auto indices_subview = Kokkos::subview(indices, team_member.team_rank(), Kokkos::ALL());
////                auto combination_subview = Kokkos::subview(combination_view, team_member.team_rank(), Kokkos::ALL());
////                auto sgraph_dist_subview = Kokkos::subview(sgraph_distance_signature, 
////                                                            team_member.team_rank(), 
////                                                            Kokkos::ALL());
////                auto sgraph_deg_subview = Kokkos::subview(sgraph_degree_signature, 
////                                                          team_member.team_rank(), 
////                                                          Kokkos::ALL());
////                auto subgraph_subview = Kokkos::subview(induced_subgraph, 
////                                                        team_member.team_rank(), 
////                                                        Kokkos::ALL(), 
////                                                        Kokkos::ALL());
////                auto visited_subview = Kokkos::subview(visited, team_member.team_rank(), Kokkos::ALL());
////                auto queue_subview = Kokkos::subview(queue, team_member.team_rank(), Kokkos::ALL());
////                auto distance_subview = Kokkos::subview(distance, team_member.team_rank(), Kokkos::ALL());
////                auto scratch_view = Kokkos::subview(scratch_block, team_member.team_rank(), Kokkos::ALL());
//          scratch_int_view scratch_view(team_member.thread_scratch(0), max_neighbors);
//          scratch_int_view indices_subview(team_member.thread_scratch(0), 5);
//          scratch_int_view combination_subview(team_member.thread_scratch(0), 5);
//          scratch_int_view sgraph_dist_subview(team_member.thread_scratch(0), orbits.distance.extent(0));
//          scratch_int_view sgraph_deg_subview(team_member.thread_scratch(0), 5);
//          scratch_int_view queue_subview(team_member.thread_scratch(0), 5);
//          scratch_int_view distance_subview(team_member.thread_scratch(0), 5);
//          scratch_subgraph_view subgraph_subview(team_member.thread_scratch(0), 5, 5);
//          scratch_bool_view visited_subview(team_member.thread_scratch(0), 5);
//                combination_from_position(scratch_view, combination_num, n_neighbors, node_count);
//                for(int j=0; j<node_count; j++) {
////                    combination_subview(j) = neighbor_subview(scratch_view(j));
//                    combination_subview(j) = neighbor_scratch(scratch_view(j));
//                }
//                combination_subview(node_count) = node;
////                kokkos_calculate_GDV(team_member, node, node_count, graph, orbits, neighbor_scratch, num_neighbors(starts(chunk_idx)), combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, metrics_sa);
////                kokkos_calculate_GDV(team_member, node, node_count, graph, orbits, neighbor_scratch, num_neighbors(starts(chunk_idx)), combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, graph_GDV);
//                kokkos_calculate_GDV(team_member, node, node_count, graph, orbits, neighbor_scratch, num_neighbors(curr_node), combination_subview, sgraph_dist_subview, sgraph_deg_subview, subgraph_subview, visited_subview, queue_subview, distance_subview, graph_GDV);
//              });
//            });
//          }
//        }
////#ifdef DEBUG
//printf("Rank %d done with chunk %d\n", rankn, chunk_idx);
////#endif
//      }
////      Kokkos::Experimental::contribute(graph_GDV, metrics_sa);
////      metrics_sa.reset();
//      Kokkos::fence();
//#ifdef ANALYSIS
//std::string region_log("Rank-");
//region_log = region_log + std::to_string(rankn) + std::string("-region-data-graph-");
//region_log = region_log + std::string(graph_name) + std::string(".log");
//std::fstream fs(region_log, std::fstream::out|std::fstream::app);
//// Hash based change tracking
//uint32_t num_changed = 0;
//#ifdef HASH_ANALYSIS
////auto test = make_layout_left(graph_GDV);
////generate_hashes(test, new_hashes, page_size);
//generate_hashes(graph_GDV, new_hashes, page_size);
//#endif
//num_changed = print_changed_blocks(fs, new_data, old_data);
//auto contiguous_regions = print_contiguous_regions(fs, new_data, old_data);
////fs << "Number of blocks: " << new_hashes.size() << std::endl;
////fs << "Changed blocks: ";
////bool flag = false;
////for(int idx = 0; idx<new_hashes.size(); idx++) {
////  if(old_hashes(idx) != new_hashes(idx))
////  {
////    num_changed += 1;
////    if(!flag) {
////      fs << "[" << idx << ",";
////      flag = true;
////    }
////  } else {
////    if(flag) {
////      fs << idx << ") ";
////      flag = false;
////    }
////  }
////}
////if(flag)
////  fs << new_hashes.size() << ")";
//////std::cout << std::endl;
//////std::cout << num_changed << "/" << new_hashes.size() << " blocks changed " << std::endl;
////fs << std::endl;
////fs << num_changed << "/" << new_hashes.size() << " blocks changed " << std::endl;
////// Find contiguous regions
////std::map<int, int> contiguous_regions;
////int largest_region = 1;
////int region_size_counter = 0;
////for(int idx=0; idx<new_hashes.size(); idx++) {
////  if(old_hashes(idx) != new_hashes(idx))
////  {
////    region_size_counter += 1;
////  }
////  else 
////  {
////    if(region_size_counter > largest_region)
////    {
////      largest_region = region_size_counter;
////    }
////    if(region_size_counter > 0) {
////      auto pos = contiguous_regions.find(region_size_counter);
////      if(pos != contiguous_regions.end())
////      {
////        pos->second = pos->second + 1;
////      }
////      else
////      {
////        contiguous_regions.insert(std::pair<int,int>(region_size_counter, 1));
////      }
////    }
////    region_size_counter = 0;
////  }
////}
////if(region_size_counter > largest_region)
////{
////  largest_region = region_size_counter;
////}
////if(region_size_counter > 0) {
////  auto pos = contiguous_regions.find(region_size_counter);
////  if(pos != contiguous_regions.end())
////  {
////    pos->second = pos->second + 1;
////  }
////  else
////  {
////    contiguous_regions.insert(std::pair<int,int>(region_size_counter, 1));
////  }
////}
////fs << "Largest region: " << largest_region << std::endl;
////fs << "Region map\n";
////for(auto itr = contiguous_regions.begin(); itr != contiguous_regions.end(); ++itr)
////{
////  fs << itr->second << " regions of size " << itr->first << std::endl;
////}
//#ifdef HASH_ANALYSIS
//Kokkos::deep_copy(old_data, new_data);
//Kokkos::deep_copy(new_data, 0);
//#endif
//#ifdef SCAN_ANALYSIS
//Kokkos::deep_copy(old_data, new_data);
//#endif
//#endif
//#ifdef DIRTY_PAGE_TRACKING
//      uintptr_t start_addr = reinterpret_cast<uintptr_t>(graph_GDV.data());
//      uintptr_t end_addr = reinterpret_cast<uintptr_t>(graph_GDV.data() + graph_GDV.span());
//      bool range_dirty = address_range_dirty(pid, start_addr, end_addr);
//      if(range_dirty) {
//        printf("Detected changes\n");
//      }
//      uint64_t num_pages = ((end_addr-start_addr)/page_size);
//      uint64_t* page_list = (uint64_t*) malloc(sizeof(uint64_t)*num_pages);
//      for(int j=0; j<num_pages; j++) {
//        page_list[j] = 0; 
//      }
//      bool dirty = get_dirty_pages(pid, start_addr, end_addr, page_list);
//      printf("List of dirty pages: ");
//      for(int j=0; j<num_pages; j++) {
//        if(page_list[j] > 0) {
//          printf("%llu, ", page_list[j]);
//        }
//      }
//      printf("\n");
//#endif
//#ifdef AUTO_CHECKPOINT
//}, filt);
//#endif
//#ifdef DETAILED_TIMERS
//      Kokkos::fence();
//      chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
//      chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2-t1);
//      printf("Done with chunk: %u, time: %f\n", chunk_idx, time_span.count());
//#endif
//      offset++;
//    }
//    if(rankn == 0) {
//      MPI_Reduce(MPI_IN_PLACE, graph_GDV.data(), graph_GDV.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
//    } else {
//      MPI_Reduce(graph_GDV.data(), graph_GDV.data(), graph_GDV.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
//    }
//  }
//  MPI_Barrier(MPI_COMM_WORLD);
//
//  vec_calc_communication_time[graph_counter - 1] = process_ends_communication - vec_calc_start;
//
//  #ifdef DEBUG
//    cout << "Finished GDV Vector Calc on Rank: " << rankn << endl;
//  #endif
//}

template<class NeighborView, class IntView, class GraphView, class BoolView, class Scatter>
KOKKOS_INLINE_FUNCTION void 
kokkos_calculate_GDV(Kokkos::TeamPolicy<>::member_type team_member,
                      int node, 
                      int node_count,
                      const matrix_type& graph, 
                      const Orbits& orbits, 
                      const NeighborView& neighbors,
                      const int num_neighbors,
                      const IntView& combination_view,
                      IntView& sgraph_distance_signature,
                      IntView& sgraph_degree_signature,
                      GraphView& induced_subgraph,
                      BoolView& visited,
                      IntView& queue,
                      IntView& distance,
                      Scatter gdvMetrics_sa
                    )
{
//  auto gdvMetrics = gdvMetrics_sa.access();
  auto combination = Kokkos::subview(combination_view, Kokkos::make_pair(0,node_count+1));
  auto subgraph_degree_signature = Kokkos::subview(sgraph_degree_signature, 
                                                    Kokkos::make_pair(0, node_count+1));
//  auto subgraph_distance_signature = Kokkos::subview(sgraph_distance_signature, 
//                                                    Kokkos::make_pair(0, orbits.distance.extent(1)));
  auto induced_sgraph = Kokkos::subview(induced_subgraph, 
                                        Kokkos::make_pair(0,node_count+1), 
                                        Kokkos::make_pair(0,node_count+1));
  EssensKokkos::kokkos_induced_subgraph(graph, combination, induced_sgraph);
  bool is_connected = EssensKokkos::is_connected(induced_sgraph, visited, queue);
  if(is_connected)
  {
    EssensKokkos::calc_degree_signature(induced_sgraph, subgraph_degree_signature);
    for(int idx=0; idx<node_count+1; idx++)
    {
      int v = idx;
      EssensKokkos::calc_distance_signature(v, induced_sgraph, sgraph_distance_signature, 
                                            visited, queue, distance);
      for(int i=orbits.start_indices(node_count+1); i<orbits.start_indices(node_count+2); i++) {
        auto orbit_deg_sig = Kokkos::subview(orbits.degree, i, Kokkos::ALL);
        bool match = EssensKokkos::compare_signatures(subgraph_degree_signature, orbit_deg_sig);
        auto orbit_dis_sig = Kokkos::subview(orbits.distance, i, Kokkos::ALL);
        match = match && EssensKokkos::compare_signatures(sgraph_distance_signature, orbit_dis_sig);
        if(match) {
//          gdvMetrics(combination(v),i) += 1;
          Kokkos::atomic_add(&gdvMetrics_sa(combination(v),i), (uint32_t)(1));
        }
      }
    }
  }
}

//template<class NeighborView, class IntView, class GraphView, class BoolView, class Scatter>
template<class IntView, class GraphView, class BoolView, class Scatter>
KOKKOS_INLINE_FUNCTION void 
kokkos_calculate_GDV2(Kokkos::TeamPolicy<>::member_type team_member,
                      int node, 
                      int node_count,
                      const matrix_type& graph, 
                      const Orbits& orbits, 
//                      const NeighborView& neighbors,
                      const int num_neighbors,
                      IntView& combination_view,
                      IntView& sgraph_distance_signature,
                      IntView& sgraph_degree_signature,
                      GraphView& induced_subgraph,
                      BoolView& visited,
                      IntView& queue,
                      IntView& distance,
                      Scatter gdvMetrics_sa
                    )
{
//  auto gdvMetrics = gdvMetrics_sa.access();
  auto combination = Kokkos::subview(combination_view, Kokkos::make_pair(0,node_count+1));
  auto subgraph_degree_signature = Kokkos::subview(sgraph_degree_signature, 
                                                    Kokkos::make_pair(0, node_count+1));
//  auto subgraph_distance_signature = Kokkos::subview(sgraph_distance_signature, 
//                                                    Kokkos::make_pair(0, orbits.distance.extent(1)));
  auto induced_sgraph = Kokkos::subview(induced_subgraph, 
                                        Kokkos::make_pair(0,node_count+1), 
                                        Kokkos::make_pair(0,node_count+1));
  EssensKokkos::kokkos_induced_subgraph(graph, combination, induced_sgraph);
  bool is_connected = EssensKokkos::is_connected(induced_sgraph, visited, queue);
  if(is_connected)
  {
    EssensKokkos::calc_degree_signature(induced_sgraph, subgraph_degree_signature);
    for(int idx=0; idx<node_count+1; idx++)
    {
      int v = idx;
      EssensKokkos::calc_distance_signature(v, induced_sgraph, sgraph_distance_signature, 
                                            visited, queue, distance);
      for(int i=orbits.start_indices(node_count+1); i<orbits.start_indices(node_count+2); i++) {
        auto orbit_deg_sig = Kokkos::subview(orbits.degree, i, Kokkos::ALL);
        bool match_deg = EssensKokkos::compare_signatures(subgraph_degree_signature, orbit_deg_sig);
        auto orbit_dis_sig = Kokkos::subview(orbits.distance, i, Kokkos::ALL);
        bool match_dis = EssensKokkos::compare_signatures(sgraph_distance_signature, orbit_dis_sig);
        if(match_deg && match_dis) {
//          gdvMetrics(combination(v),i) += 1;
          Kokkos::atomic_add(&gdvMetrics_sa(combination(v),i), static_cast<uint32_t>(1));
        }
      }
    }
  }
}


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
                    )
{
printf("Processing node %d\n", node);
//  auto gdvMetrics = gdvMetrics_sa.access();
  int num_neighbors = EssensKokkos::find_neighbours(node, graph, 4, neighbor_buff);
  auto neighbors = Kokkos::subview(neighbor_buff, Kokkos::make_pair(0, num_neighbors));
printf("Found neighbors\n");
  for (int node_count = 1; node_count < 5; node_count++)
  {
    CombinationGenerator generator(num_neighbors, node_count, indices);
    while(!generator.done) 
    {
      auto combination = Kokkos::subview(combination_view, Kokkos::make_pair(0,node_count+1));
      generator.kokkos_get_combination(indices, num_neighbors, neighbors, combination);
      auto subgraph_degree_signature = Kokkos::subview(sgraph_degree_signature, 
                                                        Kokkos::make_pair(0, node_count+1));
      auto subgraph_distance_signature = Kokkos::subview(sgraph_distance_signature, 
                                                        Kokkos::make_pair(0, static_cast<int>(orbits.distance.extent(1))));
      combination(node_count) = node;
printf("Generated combination for node %d\n", node);
      auto induced_sgraph = Kokkos::subview(induced_subgraph, 
                                            Kokkos::make_pair(0,node_count+1), 
                                            Kokkos::make_pair(0,node_count+1));
      EssensKokkos::kokkos_induced_subgraph(graph, combination, induced_sgraph);
printf("Created induced subgraph for node %d\n", node);
      bool is_connected = EssensKokkos::is_connected(induced_sgraph, visited, queue);
printf("Checked if graph is connected for node %d\n", node);
printf("Graph for %d is connected\n", node);
      if(is_connected)
      {
        EssensKokkos::calc_degree_signature(induced_sgraph, subgraph_degree_signature);
        for(int idx=0; idx<node_count+1; idx++)
        {
          int v = idx;
          EssensKokkos::calc_distance_signature(v, induced_sgraph, subgraph_distance_signature, 
                                                visited, queue, distance);
          for(int i=orbits.start_indices(node_count+1); i<orbits.start_indices(node_count+2); i++) 
          {
            auto orbit_deg_sig = Kokkos::subview(orbits.degree, i, Kokkos::ALL);
            bool match = EssensKokkos::compare_signatures(subgraph_degree_signature, orbit_deg_sig);
            auto orbit_dis_sig = Kokkos::subview(orbits.distance, i, Kokkos::ALL);
            match = match && EssensKokkos::compare_signatures(subgraph_distance_signature, orbit_dis_sig);
            if(match) 
            {
//              gdvMetrics(combination(v),i) += 1;
              Kokkos::atomic_add(&gdvMetrics_sa(combination(v),i), static_cast<uint32_t>(1));
            }
          }
        }
      }
      generator.kokkos_next(indices);
    }
printf("Completed subgraphs of size %d for node %d\n", node_count, node);
  }
}

