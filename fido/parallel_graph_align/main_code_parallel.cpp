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
#include <time.h>
#include <stdlib.h>
#include <ctime>
#include <mpi.h>
#include <fstream>
#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include <Kokkos_ScatterView.hpp>

#define GDV_LENGTH 22
#define CHUNK_SIZE 1
//#define DEBUG
//#define INSTRUMENT

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
KOKKOS_INLINE_FUNCTION void kokkos_calculate_GDV(int node, const matrix_type& graph, const Orbits& orbits, GDVs& gdvMetrics);
KOKKOS_INLINE_FUNCTION void kokkos_calculate_GDV(int node, const matrix_type& graph, const Orbits& orbits, Kokkos::Experimental::ScatterView<int**> gdvMetrics);

struct node_id_order {
  inline bool operator() (const GDVMetric& gdv_met_1, const GDVMetric& gdv_met_2) {
    return (gdv_met_1.node < gdv_met_2.node);
  }
};

// Define variables for keeping track of time for load imbalancing tests.
double total_time_taken;
double vec_calc_communication_time[2] = {};
double* vec_calc_computation_time_X;
double* vec_calc_computation_time_Y;
int* vec_calc_proc_assign_X;
int* vec_calc_proc_assign_Y;
//double vec_calc_prior_gather;
//double vec_calc_post_gather;
int vec_calc_avg_node_deg;

using namespace std;

int main(int argc, char *argv[]) {

  MPI_Init(&argc,&argv);
//  int provided;
//  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided); 
  Kokkos::initialize(argc, argv);
  {

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

//  ifstream the_file3 ( argv[4] + "time_results.txt" );
//  if (!the_file3.is_open() ) {
//    cout << "INPUT ERROR:: Could not open the time recording file\n";
//  }

//  vector<OrbitMetric> orbits;
//  readin_orbits(&the_file2,&orbits);
  Orbits k_orbits;
  kokkos_readin_orbits(&the_file2, k_orbits);

  A_Network X;
  A_Network Y;
  readin_network(&X,argv[1],0,-1);
  readin_network(&Y,argv[2],0,-1);
  GDV_functions test_gdvf;

  matrix_type graphx, graphy;
  readin_graph(&the_file0, graphx);
  readin_graph(&the_file1, graphy);

for(int i=0; i<X.size(); i++) {
  auto row = graphx.rowConst(i);
//cout << "A_Network: " << X[i].ListW.size() << " Graph: " << row.length << " row: " << X[i].Row << " " << i << endl;
  for(int j=0; j<X[i].ListW.size(); j++) {
    int xv = X[i].ListW[j].first;
    int gxv = row.colidx(j);
    if(xv != gxv) {
      cout << "Graphs don't match\n";
      cout << xv << " vs " << gxv << endl;
      break;
    }
  }
}
for(int i=0; i<Y.size(); i++) {
  auto row = graphy.rowConst(i);
//cout << "A_Network: " << X[i].ListW.size() << " Graph: " << row.length << " row: " << X[i].Row << " " << i << endl;
  for(int j=0; j<Y[i].ListW.size(); j++) {
    int xv = Y[i].ListW[j].first;
    int gxv = row.colidx(j);
    if(xv != gxv) {
      cout << "Graphs don't match\n";
      cout << xv << " vs " << gxv << endl;
      break;
    }
  }
}

  //clock_t out_tStart = clock();

  int numtasks, rank, dest, source, rc, count, tag=0;
  MPI_Status Stat;   // required variable for receive routines                                                                                                                                                                          

  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

  // Perform Similarity Calculations
//  Similarity_Metric_calculation_for_two_graphs(X,Y,orbits, graph_name1, graph_name2);
  kokkos_Similarity_Metric_calculation_for_two_graphs(graphx, graphy, k_orbits, graph_name1, graph_name2);

  #ifdef DEBUG
    cout << "Finished Similarity Metric Calculation on Rank: " << rank << endl;
  #endif

  // Set up for and Perform Runtime Management and Gather
  double* time_buff = NULL;
  int num_times = 3;
  double send_times[num_times];
//  if (rank == 0) {
//    vec_calc_computation_time = 0;
//  }
//  cout << "Runtime on rank " << rank << " for graph 1 = " << vec_calc_communication_time[0] << endl;
//  cout << "Runtime on rank " << rank << " for graph 2 = " << vec_calc_communication_time[1] << endl;
  send_times[0] = total_time_taken;
  send_times[1] = vec_calc_communication_time[0];
  send_times[2] = vec_calc_communication_time[1];
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
      myfile << month << "/" << day << "/" << year << "\t" << numtasks << "\t" << graph_name1 << "\t" << graph_name2 << "\t" << graphx.numRows() << "\t\t" << total_time_taken << " \n";
      myfile.close();
    }
    else { 
      cout << "Out File Did Not Open";
    }

    // Print out rank specific runtime data
    string time_file = "runtime_data/runtimes_rec_over.txt";
    myfile.open(time_file, ofstream::trunc);
    if (!myfile.is_open()) {
      cout << "INPUT ERROR:: Could not open the local time recording file\n";
    }
    myfile << "Time Taken in Similarity Metric Calculation = " << " \n";
    myfile << "Rank\tGraph 1\tGraph 2\tTotal\n";
    for (int i = 0; i < numtasks; i++) {
      myfile << i << " " << time_buff[num_times*i+1] << " " << time_buff[num_times*i+2] << " " << time_buff[num_times*i] << " \n";
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
//cout << "Number of orbits: " << num_orbits << ", number of vertices: " << graph1.numRows() << endl;
  GDVs graph1_GDV("GDV for graph1", graph1.numRows(), num_orbits);
  GDVs graph2_GDV("GDV for graph2", graph2.numRows(), num_orbits);

  int rankm, numtasksm;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankm);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasksm);

  int graph_counter = 1;
  kokkos_GDV_vector_calculation(graph1, graph1_GDV, orbits, "graph1", graph_counter); 

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
  MPI_Barrier(MPI_COMM_WORLD);

  int m = graph1_GDV.extent(0);
  int n = graph2_GDV.extent(0);
  Kokkos::View<double**> sim_mat("Similarity matrix", m, n);
//  Kokkos::MDRangePolicy<Kokkos::Rank<2>> mdrange({0,0}, {m,n});
//  Kokkos::parallel_for("Compute sim_mat", mdrange, KOKKOS_LAMBDA(const int i, const int j) {
if(rankm == 0) {
  for(int i=0; i<m; i++) {
    for(int j=0; j<n; j++) {
      auto g1_subview = Kokkos::subview(graph1_GDV, i, Kokkos::ALL);
      auto g2_subview = Kokkos::subview(graph2_GDV, j, Kokkos::ALL);
      sim_mat(i,j) = kokkos_GDV_distance_calculation(g1_subview, g2_subview);
    }
  }
}
//  });
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
        myfile<<" { "<<sim_mat(i,j)<<" } ";
      }
      myfile<<"||"<<endl;
    }
    myfile.close();
    
    
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
  int sum = 0;
  for(int i=0; i<gdvm.extent(0); i++) {
    sum += gdvm[i]*gdvm[i];
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
  int graph_size = graph.numRows();
cout << "Starting GDV vector calculation\n";
//  int graph_size = graph.size();
//  graph_GDV->clear();
//
//  for(int i=0; i<graph.size(); i++) {
//    vector<int> gdv_vec(GDV_LENGTH, 0);
//    graph_GDV->push_back(GDVMetric(graph[i].Row, gdv_vec));
//  }
  //graph_counter += 1;
  //cout << "Graph counter on rank " << rankn << " = " << graph_counter << endl;

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

      process_ends_communication = MPI_Wtime();
      //vec_calc_prior_gather = MPI_Wtime() - vec_calc_start + vec_calc_prior_gather;

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
        process_ends_communication = MPI_Wtime();
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
#ifdef INSTRUMENT
printf("Combination | # of updated regions | avg size of region | (# of identical regions, total size)\n");
#endif
      GDVs metrics("GDVs", graph.numRows(), GDV_LENGTH);
//cout << "created GDV\n";
      Kokkos::Experimental::ScatterView<int**> metrics_sa(metrics);
//Kokkos::parallel_for("Calculate GDV", Kokkos::RangePolicy<>(node_name,end), KOKKOS_LAMBDA(const int node) {
//cout << "Entering kernel\n";
      for(int node=node_name; node<end; node++) {
        vec_calc_computation_start = MPI_Wtime();
//cout << "Starting gdv calculation\n";
//        Calculate_GDV(graph[node].Row, graph, orbits, metrics);
        kokkos_calculate_GDV(node, graph, orbits, metrics_sa);
//        kokkos_calculate_GDV(node, graph, orbits, metrics);

        if (graph_counter == 1) {
          vec_calc_computation_time_X[node] = MPI_Wtime() - vec_calc_computation_start;
          vec_calc_proc_assign_X[node] = rankn;
        }
        if (graph_counter == 2) {
          vec_calc_computation_time_Y[node] = MPI_Wtime() - vec_calc_computation_start;
          vec_calc_proc_assign_Y[node] = rankn;
        }
        #ifdef DEBUG
//          cout << "Sending GDV for node " << graph[node_name].Row << " from rank " << rankn << endl;
          cout << "Sending GDV for node " << node_name << " from rank " << rankn << endl;
        #endif
      }
//});
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
  }

  //vec_calc_post_gather = MPI_Wtime() - vec_calc_start + vec_calc_post_gather;
  vec_calc_communication_time[graph_counter - 1] = process_ends_communication - vec_calc_start;
//  cout << "Communication time on process " << rankn << " for graph " << graph_counter << " = " << vec_calc_communication_time[graph_counter - 1] << endl;
  //vec_calc_computation_time = vec_calc_computation_end - vec_calc_computation_start + vec_calc_computation_time;

  #ifdef DEBUG
    cout << "Finished GDV Vector Calc on Rank: " << rankn << endl;
  #endif

}

KOKKOS_INLINE_FUNCTION void 
//void kokkos_calculate_GDV(int node, const matrix_type& graph, const Orbits& orbits, GDVs& gdvMetrics)
kokkos_calculate_GDV(int node, 
                      const matrix_type& graph, 
                      const Orbits& orbits, 
                      Kokkos::Experimental::ScatterView<int**> gdvMetrics_sa)
{
  auto gdvMetrics = gdvMetrics_sa.access();
  Kokkos::View<int*> neighbors;
  EssensKokkos::find_neighbours(node, graph, 4, neighbors);
//cout << "Allocated and found neighbors of " << node << ": " << neighbors.size() << " total neighbors\n";
//for(int i=0; i<neighbors.size(); i++) {
//  cout << neighbors(i) << " ";
//}
//cout << endl;
  Kokkos::View<int*> neighbor_set("Set", neighbors.size());
  Kokkos::deep_copy(neighbor_set, neighbors);
  for (int node_count = 1; node_count < 5; node_count++)
  {
//cout << "Starting " << node_count << " node subgraphs\n";
    CombinationGenerator generator(neighbors.size(), node_count);
    Kokkos::View<int*> combination("Combination", node_count+1);
#ifdef INSTRUMENT
Kokkos::View<bool**> gdv_updated("Updated entries", graph.numRows(), orbits.distance.extent(0));
Kokkos::View<bool**> gdv_updated_prev("Previously updated entries", graph.numRows(), orbits.distance.extent(0));
Kokkos::deep_copy(gdv_updated, false);
Kokkos::deep_copy(gdv_updated_prev, false);
#endif
    while(!generator.done) {
      generator.kokkos_get_combination(neighbors, combination);
      matrix_type induced_sgraph;
      Kokkos::View<int*> subgraph_degree_signature("Subgraph degree signature", node_count+1);
      Kokkos::View<int*> subgraph_distance_signature("Subgraph distance signature", orbits.distance.extent(1));
      combination(node_count) = node;
      EssensKokkos::kokkos_inducedSubgraph(graph, combination, induced_sgraph);
      bool is_connected = EssensKokkos::isConnected(induced_sgraph);
      if(is_connected)
      {
//if(node_count == 1) {
//cout << "Found valid combination: ";
//for(int i=0; i<combination.size(); i++) {
//  cout << combination(i) << " ";
//}
//cout << endl;
//}
        EssensKokkos::degree_signature(induced_sgraph, subgraph_degree_signature);
        Kokkos::sort(subgraph_degree_signature);
//if(node_count == 1) {
//cout << "Degree signature: ";
//for(int i=0; i<subgraph_degree_signature.size(); i++) {
//  cout << subgraph_degree_signature(i) << " ";
//}
//cout << endl;
//}
        for(int idx=0; idx<induced_sgraph.numRows(); idx++)
        {
          int v = idx;
          EssensKokkos::distance_signature(v, induced_sgraph, subgraph_distance_signature);
//if(node_count == 1) {
//cout << "Distance signature: ";
//for(int i=0; i<subgraph_distance_signature.size(); i++) {
//  cout << subgraph_distance_signature(i) << " ";
//}
//cout << endl;
//}
          for(int i=orbits.start_indices(node_count+1); i<orbits.start_indices(node_count+2); i++) {
            bool match = true;
            for(int j=0; j<node_count+1; j++) {
              if(subgraph_degree_signature(j) != orbits.degree(i,j)) {
                match = false;
                break;
              }
            }
            if(!match) {
              continue;
            }
            for(int j=0; j<orbits.distance.extent(1); j++) {
              if(subgraph_distance_signature(j) != orbits.distance(i,j)) {
                match = false;
                break;
              }
            }
            if(match) {
//cout << "Updating gdv\n";
              gdvMetrics(combination(v),i) += 1;
#ifdef INSTRUMENT
gdv_updated(combination(v), i) = true;
#endif
            }
          }
        }
#ifdef INSTRUMENT
cout << "(";
for(int i=0; i<combination.size()-1; i++) {
  printf("%4d, ", combination(i));
}
printf("%4d)", combination(combination.size()-1));
for(int i=0; i<5-combination.size(); i++) {
  printf("      ");
}
cout << "\t\t";
bool contiguous = false;
int gdv_start = 0;
int gdv_end = 0;
int region_start = 0;
int region_end = 0;
int region_size = 0;
int num_regions = 0;
for(int i=0; i<gdv_updated.extent(0); i++) {
  for(int j=0; j<gdv_updated.extent(1); j++) {
    if(gdv_updated(i,j) && !contiguous) {
      contiguous = true;
      gdv_start = i;
      region_start = j;
      region_size += 1;
    } else if(gdv_updated(i,j) && contiguous) {
      region_size += 1;
    } else if(!gdv_updated(i,j) && contiguous) {
      contiguous = false;
      gdv_end = i;
      region_end = j;
      num_regions += 1;
    }
  }
}
cout << "(" << num_regions << ")\t";
if(num_regions == 0) {
  printf("(%7f)\t", 0.0);
} else {
  printf("(%7f)\t", static_cast<float>(region_size)/static_cast<float>(num_regions));
}
region_size = 0;
num_regions = 0;
contiguous = false;
for(int i=0; i<gdv_updated.extent(0); i++) {
  for(int j=0; j<gdv_updated.extent(1); j++) {
    if(gdv_updated(i,j) && gdv_updated_prev(i,j) && !contiguous) {
      contiguous = true;
      region_size += 1;
    } else if(gdv_updated(i,j) && gdv_updated_prev(i,j) && contiguous) {
      region_size += 1;
    } else if((!gdv_updated(i,j) || !gdv_updated_prev(i,j)) && contiguous) {
      contiguous = false;
      num_regions += 1;
    }
  }
}
cout << "(" << num_regions << "," << region_size << ")";
Kokkos::deep_copy(gdv_updated_prev, gdv_updated);
cout << endl;
#endif
      }
      generator.kokkos_next();
    }
  }
//cout << "Inside kokkos_calculate_GDV: GDV: Root node " << node << endl;
//for(int i=0; i<gdvMetrics.extent(0); i++) {
//  cout << "Node " << i << ": ";
//  for(int j=0; j<gdvMetrics.extent(1); j++) {
//    cout << gdvMetrics(i,j) << " ";
//  }
//  cout << endl;
//}
//cout << endl;
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

  int lastrow;
  vector<int> rows, cols;
  vector<float> vals;
  int k = 0;
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
if(w != 0.0) {
    edge1.push_back(1);
} else {
    edge1.push_back(0);
}
//    edge1.push_back(w);
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
//      edge2.push_back(w);
      edge_list.push_back(edge2);
    }
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

