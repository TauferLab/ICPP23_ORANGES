// Required directives are to be added from ESSENS to use the functionality 

#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
#include "GDV_functions.hpp"

using namespace std;

void GDV_functions::find_combination ( int set[], int size, vector<vector<int>>& output) {
  
  return;

};

void GDV_functions::inducedSubGraph(A_Network network, vector<int> nodes, A_Network& output) {
    
  return;

}

void GDV_functions::isConnected(A_Network network, bool& isConnected) {

  return;

}

void GDV_functions::find_neighbours(int node,A_Network network,int distance,vector<int> &neighbours)
//GDV_functions::find_neighbours(int node, vector<ADJ_Bundle> network, int distance, vector<int> &neighbours) 
{
  
  get_unq_neighbors(node, network, distance, &neighbours);
  return;

}

void GDV_functions::test_func() {
}
