// Required directives are to be added from ESSENS to use the functionality 

#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
#include "GDV_functions.hpp"
#include "class_definitions.hpp"
GDVMetric Calculate_GDV(int ,A_Network );
using namespace std;

int main(int argc, char *argv[]) {
       
  clock_t q, q1, q2,t;
 
  /* Accepts file as input. 
     Input should be in the format of ( node1 node2 weight ). 
     By default, weight should be 1
     Should be sorted. */ 
  ifstream the_file1 ( argv[1] ); 
  if (!the_file1.is_open() ) { 
    cout<<"INPUT ERROR:: Could not open the input file\n";
  }
  
  /* Read in Network: Reads the file and converts it into a network of type A_Network*/
  A_Network X;
  readin_network(&X,argv[1],1,-1);  
 

  // Objects for testing GDV induced subgraph function
  // A_Network subgraph;
  // vector<int> subgraph_nodes;
  // subgraph_nodes.push_back(0);
  // subgraph_nodes.push_back(1);
  // subgraph_nodes.push_back(2);
  // subgraph_nodes.push_back(6);
  // gdvf.inducedSubgraph(X, subgraph_nodes, subgraph);
  //gdvf.print_disconnected_network(subgraph);

  // Objects for testing connectedness function
  // bool is_connected = false;
  // gdvf.isConnected(subgraph, is_connected);
  //cout << is_connected << endl;

    
  // for (int i=0;i<X.size(); i++)
  // {
  //   // Calculate_GDV(i,X);
  // }
  vector<int> GDV;
  int node;
  GDVMetric gdv(node,GDV);
  gdv = Calculate_GDV(2,X);
  return 0;
}

GDVMetric Calculate_GDV(int node,A_Network Graph)
{
    GDV_functions gdvf;
    printf("calculating GDV for node %d\n",node);
    vector<int> neighbours;
    gdvf.find_neighbours(node,Graph,2,&neighbours);
    int set[neighbours.size()]; 
    std::copy( neighbours.begin(), neighbours.end(), set );
    int numElements = *(&set + 1) - set;
    vector<vector<int>> combinationsList;
    // // for info 
    // for (int i = numElements - 1; i >= 0; i--) 
    // {
    //   cout << "element "<<i<< ":"<<set[i] << endl;
    // }
    for (int node_count = 0; node_count < 5; node_count++)
    {
      gdvf.find_combinations(set, numElements,node_count,&combinationsList);
      for (vector<int> combination : combinationsList)
      {
        A_Network induced_sgraph;
        bool is_connected = false;
        combination.push_back(node);
        print_vector(combination);
        gdvf.inducedSubgraph(Graph, combination, induced_sgraph);
        cout<<"sub graph for combination "<<endl;
        print_network(induced_sgraph);
        gdvf.print_disconnected_network(induced_sgraph);
        gdvf.isConnected(induced_sgraph, is_connected);
        if(is_connected)
        {
          cout<<"into connected"<<endl;
        }
      }
    }

}



