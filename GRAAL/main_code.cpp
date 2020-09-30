// Required directives are to be added from ESSENS to use the functionality 

//#include "structure_defs.hpp"
//#include "input_to_network.hpp"
//#include "printout_others.hpp"
//#include "printout_network.hpp"
//#include "ADJ/find_Xneighbors.hpp"
#include "GDV_functions.hpp"

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
  GDV_functions gdvf;

  // Objects for testing GDV induced subgraph function
  A_Network subgraph;
  vector<int> subgraph_nodes;
  subgraph_nodes.push_back(1);
  subgraph_nodes.push_back(2);
  subgraph_nodes.push_back(3);
  //gdvf.inducedSubgraph(X, subgraph_nodes, subgraph);
  //print_network(subgraph);

  // Objects for testing connectedness function
  bool is_connected = false;
  //gdvf.isConnected(subgraph, &is_connected);

  // looping through nodes 
  // Loop through all the vertices. currently considered integers starting from 0 to n and none are missing.
  for (int i=0;i<X.size(); i++ ) {
    printf("%d\n",i);
    vector<int> neighbours;
    gdvf.find_neighbours(i,X,2,&neighbours);
    print_vector(neighbours);
    
  }

  print_network(X);

  return 0;

}



