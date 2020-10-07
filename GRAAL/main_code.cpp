// Required directives are to be added from ESSENS to use the functionality 

#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
#include "GDV_functions.hpp"
#include "class_definitions.hpp"

GDVMetric Calculate_GDV(int ,A_Network );
void readin_orbits(  ifstream* ,vector<OrbitMetric>* );
void convert_string_vector_int(string* , vector<int>* ,string );
using namespace std;

int main(int argc, char *argv[]) {
       
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

  //  ifstream the_file2 ( argv[2] ); 
  // if (!the_file2.is_open() ) { 
  //   cout<<"INPUT ERROR:: Could not open the orbit file\n";
  // }
  
  // vector<OrbitMetric> orbits;
  // readin_orbits(&the_file2,&orbits);

  /* Read in Network: Reads the file and converts it into a network of type A_Network*/
  A_Network X;
  readin_network(&X,argv[1],1,-1);  
  GDV_functions test_gdvf;

  // Objects for testing GDV induced subgraph function
  // A_Network subgraph;
  // vector<int> subgraph_nodes;
  // subgraph_nodes.push_back(0);
  // subgraph_nodes.push_back(1);
  // subgraph_nodes.push_back(2);
  // subgraph_nodes.push_back(6);
  // gdvf.inducedSubgraph(X, subgraph_nodes, subgraph);
  // // gdvf.print_disconnected_network(subgraph);
  // cout<<"subgraph for 0,1,2,6"<<endl;
  // print_network(subgraph);

  // Objects for testing connectedness function
  // bool is_connected = false;
  // gdvf.isConnected(subgraph, is_connected);
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

  vector<int> GDV;
  int node;
  //GDVMetric gdv(node,GDV);
  // Calculate_GDV(1,X);
  return 0;
}

GDVMetric Calculate_GDV(int node,A_Network Graph)
{
    GDV_functions gdvf;
    printf("calculating GDV for node %d\n",node);
    vector<int> neighbours;
    gdvf.find_neighbours(node,Graph,4,&neighbours);
    print_vector(neighbours);
    int set[neighbours.size()]; 
    std::copy( neighbours.begin(), neighbours.end(), set );
    int numElements = *(&set + 1) - set;
    
    // // for info 
    // for (int i = numElements - 1; i >= 0; i--) 
    // {
    //   cout << "element "<<i<< ":"<<set[i] << endl;
    // }
    for (int node_count = 1; node_count < 5; node_count++)
    {
      vector<vector<int>> combinationsList;
      gdvf.find_combinations(set, numElements,node_count,&combinationsList);
      cout<<"Node count is "<<node_count<<endl;
      cout<<"total combinations are : "<<combinationsList.size()<<endl;
      for (vector<int> combination : combinationsList)
      {
        A_Network induced_sgraph;
        bool is_connected = false;
        combination.push_back(node);
        for (int element:combination )
        {
          cout<<" | "<<element;
        }
        cout<<endl;
        cout<<"----------------"<<endl;
        // print_vector(combination);
        // gdvf.inducedSubgraph(Graph, combination, induced_sgraph);
        // cout<<"sub graph for combination "<<endl;
        // print_network(induced_sgraph);
        // // gdvf.print_disconnected_network(induced_sgraph);
        // gdvf.isConnected(induced_sgraph, is_connected);
        // if(is_connected)
        // {
        //   cout<<"into connected"<<endl;
        // }
      }
    }

}

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
   
    
    while ((pos = s.find(signature_delimiter)) != std::string::npos) 
    {
      vector<int> segment; 
      string token;
      token = s.substr(0, pos);
      token.erase(remove(token.begin(), token.end(), '['), token.end());
      token.erase(remove(token.begin(), token.end(), ']'), token.end());
      cout<<"token is : "<<token<<endl;
      convert_string_vector_int(&token,&segment,internal_delimiter);
      s.erase(0, pos + signature_delimiter.length());
      vector_line.push_back(segment);
    }

    OrbitMetric orbMetric(vector_line[0][0],vector_line[1],vector_line[2]);
    cout<<"orbit "<<vector_line[0][0]<<endl;
    cout<<"degree signature is "<<endl;
    print_vector(vector_line[1]);
    cout<<"distance signature is "<<endl;
    print_vector(vector_line[2]);
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



