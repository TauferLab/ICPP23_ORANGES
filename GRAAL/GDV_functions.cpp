// Required directives are to be added from ESSENS to use the functionality 

#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
#include "GDV_functions.hpp"

using namespace std;

// Constructs lists of combinations of set and stores the combinations in output.
void GDV_functions::find_combination ( int set[], int size, vector<vector<int>>& output) {
  
  return;
};

// Constructs a graph out of nodes from network and stores the new graph in output.
void GDV_functions::inducedSubGraph(A_Network network, vector<int> nodes, A_Network& output) {
    
  return;
}

// Evaluates whether network is a connected graph and stores the result in isConnected.
void GDV_functions::isConnected(A_Network network, bool& isConnected) {

  return;
}

// Evaluates the neighbors of node 'node' in graph 'network' up to the given distance.  List of neighbors are stored in neighbours.
void GDV_functions::find_neighbours(int node,A_Network network,int distance,vector<int> &neighbours) {
  
  get_unq_neighbors(node, network, distance, &neighbours);
  return;
}
