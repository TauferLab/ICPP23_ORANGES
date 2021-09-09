// Required directives are to be added from ESSENS to use the functionality 

#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
#include "headers/GDV_functions.hpp"
#include "headers/class_definitions.hpp"
#include "headers/print_disconnected_graph.hpp"
#include <time.h>
#include <math.h>
#include <omp.h>

void Calculate_GDV(int ,A_Network ,vector<OrbitMetric>&, GDVMetric&);
void readin_orbits(  ifstream* ,vector<OrbitMetric>* );
void convert_string_vector_int(string* , vector<int>* ,string );
void GDV_vector_calculation(A_Network graph,vector<GDVMetric>* graph_GDV,  vector<OrbitMetric> orbits,int p);
void metric_formula(GDVMetric gdvm, double* gdv_score);
double GDV_distance_calculation(GDVMetric gdvm1, GDVMetric gdvm2);
void Similarity_Metric_calculation_for_two_graphs(A_Network graph1, A_Network graph2,vector<OrbitMetric> orbits,int p);
using namespace std;

class graph_combination
{
   public:
   int parent;
   int child;
};

class subgraph
{
    public:
    int root;
    vector<graph_combination> list;
    vector<int> list_of_nodes;
    vector<int> end_nodes;
    A_Network adj_list;
    bool is_subgraph;
};
void generate_combinations(vector<graph_combination> &list_to_generate_combinations ,vector<vector<graph_combination>> &combination_list)
{
    long int size_of_list = list_to_generate_combinations.size();
    long int iterations = pow(2,size_of_list) - 1;

    for(long int i = 1; i <= iterations; i++)
    {
        vector<graph_combination> combination_list_temp;
        int track_bit = size_of_list - 1;
        long int j = i; 
        while(j > 0)
        {
            if(j&1 == 1)
            {
                combination_list_temp.push_back(list_to_generate_combinations[track_bit]);
            }
            track_bit--;
            j = j>>1;
        }
        combination_list.push_back(combination_list_temp);
    }

}
void subgraph_enumeration(A_Network graph, vector<subgraph> &subgraph_combinations)
{
    for(int i = 0; i < graph.size(); i++)
    {
        ADJ_Bundle node = graph[i];
        //vector<subgraph> sub_graph_combinations;
        for(int j=0; j<node.ListW.size(); j++)
        {
            subgraph temp_subgraph;
            // making root
            temp_subgraph.root = i;
            // making parent-child combination
            graph_combination temp_list;
            temp_list.parent = i;
            temp_list.child = node.ListW[j].first;
            temp_subgraph.list.push_back(temp_list);
            // making end nodes
            temp_subgraph.end_nodes.push_back(node.ListW[j].first);
            temp_subgraph.list_of_nodes.push_back(i);
            temp_subgraph.list_of_nodes.push_back(node.ListW[j].first);
            making adjacency list
            ADJ_Bundle item;
            item.Row = i;
            int_double temp_ListW;
            temp_ListW.first = node.ListW[j].first;
            temp_ListW.second = 1;
            item.ListW.push_back(temp_ListW);
            temp_subgraph.adj_list.push_back(item);
            if(node.ListW[j].first > i)
            {
                temp_subgraph.is_subgraph = 1;
            }
            else
            {
                temp_subgraph.is_subgraph = 0;
            }
            subgraph_combinations.push_back(temp_subgraph);
        }

    }

    for(subgraph element1 : subgraph_combinations)
    {
        vector<vector<graph_combination>> combination_list;
        vector<graph_combination> temp_combination;
        for(int k= 0; k< element1.end_nodes.size(); k++)
        {
            int end_node_item;
            end_node_item = element1.end_nodes[k];
            ADJ_Bundle node_temp = graph[element1.end_nodes[k]];
            for(int p = 0; p<node_temp.ListW.size(); p++)
            {
              graph_combination temp_comb;
              temp_comb.parent = element1.end_nodes[k];
              temp_comb.child = node_temp.ListW[p].first;
              temp_combination.push_back(temp_comb);
            }
        }
        generate_combinations(temp_combination, combination_list);
        for(vector<graph_combination> element : combination_list)
       {
        for(graph_combination element_1 : element)
        {
            cout << element_1.parent;
            cout <<"||";
            cout << element_1.child;
            cout <<" ";
        }
        cout <<endl;
      }

    }

} 

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
    cout<<"INPUT ERROR:: Could not open the second graph input file\n";
  }

   ifstream the_file3 ( argv[3] ); 
  if (!the_file3.is_open() ) { 
    cout<<"INPUT ERROR:: Could not open the orbit file\n";
  }
      int p = atoi(argv[4]);
  vector<OrbitMetric> orbits;
  readin_orbits(&the_file3,&orbits);
  // Objects for testing orbit creation 
  // print_vector(orbits[1].orbitDegree);
  // vector<OrbitMetric> filter_o = gdvf.orbit_filter(&orbits,3);
  /* Read in Network: Reads the file and converts it into a network of type A_Network*/
  A_Network X;
  readin_network(&X,argv[1],0,-1);  


  A_Network Y;
  readin_network(&Y,argv[2],0,-1);  
  print_network(X);
  print_network(Y);
  vector<subgraph> subgraph_combinations;
  subgraph_enumeration(X,subgraph_combinations);
  cout<<"*************"<<endl;
  cout<<subgraph_combinations.size();
    /*for(vector<subgraph> element : subgraph_combinations)
    {
        print_network(element.adj_list);
    }*/
  //print_network(subgraph_combinations[0].adj_list);
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
  //cout << is_connected << endl;

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

 //Similarity_Metric_calculation_for_two_graphs(X,Y,orbits,p);
 printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  return 0;
}

void Similarity_Metric_calculation_for_two_graphs(A_Network graph1, A_Network graph2,vector<OrbitMetric> orbits ,int p)
{
  vector<GDVMetric> graph1_GDV;
  vector<GDVMetric> graph2_GDV;
 double start = omp_get_wtime();
  GDV_vector_calculation(graph1, &graph1_GDV, orbits,p); 
  GDV_vector_calculation(graph2, &graph2_GDV,orbits,p); 
  
  // vector<vector<int>> similarity_matrix( graph1_GDV.size() , vector<int> (graph2_GDV.size(), 0));
  int m = (int)graph1_GDV.size();
  int n = (int)graph2_GDV.size();

  double sim_mat[m][n];
  for (GDVMetric gdvm1: graph1_GDV)
  {
    for(GDVMetric gdvm2: graph2_GDV)
    {
      sim_mat[gdvm1.node][gdvm2.node]= GDV_distance_calculation(gdvm1,gdvm2);
    }
  }
double end = omp_get_wtime();
cout<<"Omp time taken is " <<end - start<<endl;
  ofstream myfile;
  string filename; 
  filename = "out_similarity_matrix.txt";
  myfile.open(filename);

  for(int i=1; i<m;i++)
  {
    for(int j=1;j<n;j++)
    {
      myfile<<" { "<<sim_mat[i][j]<<" } ";
    }
    myfile<<"||"<<endl;
  }
  myfile.close();

}

double GDV_distance_calculation(GDVMetric gdvm1, GDVMetric gdvm2)
{   
  double gdv1_score;
  double gdv2_score;
  metric_formula(gdvm1,&gdv1_score);
  metric_formula(gdvm2,&gdv2_score);
  
  double similarity_score = abs(gdv1_score - gdv2_score);

  return similarity_score;
}

void metric_formula(GDVMetric gdvm, double* gdv_score)
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

void GDV_vector_calculation(A_Network graph,vector<GDVMetric>* graph_GDV,  vector<OrbitMetric> orbits,int p)
{
	
omp_set_num_threads(p);
#pragma omp parallel for num_threads(p) schedule(static)
for(int i=0;i<graph.size();i++)
{
	ADJ_Bundle node = graph[i];
	vector<int> GDV_1;
     GDVMetric gdvMetric(node.Row,GDV_1);
     Calculate_GDV(node.Row,graph,orbits,gdvMetric);
     // cout<<"gdv for node "<<node.Row<<endl;
     graph_GDV->push_back(gdvMetric);
}
// for (ADJ_Bundle node:graph)
  // {
    // vector<int> GDV_1;
    // GDVMetric gdvMetric(node.Row,GDV_1);
    // Calculate_GDV(node.Row,graph,orbits,gdvMetric);
    // // cout<<"gdv for node "<<node.Row<<endl;
    // graph_GDV->push_back(gdvMetric);
  // }
}
void Calculate_GDV(int node,A_Network Graph,vector<OrbitMetric> &orbits, GDVMetric &gdvMetric)
{
    GDV_functions gdvf;
    vector<int> gdv(orbits.size(),0);
    // printf("calculating GDV for node %d\n",node);
    vector<int> neighbours;
    gdvf.find_neighbours(node,Graph,4,&neighbours);
    // print_vector(neighbours);
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



