// Required directives are to be added from ESSENS to use the functionality 

#ifndef GDV_FUNCTIONS_HPP
#define GDV_FUNCTIONS_HPP

#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
using namespace std;

class GDV_functions{

public:
    // Add the implementation over here for now. 
    void find_combination ( int set[], int size, vector<vector<int>> *output){return;}
    void inducedSubGraph(A_Network network, vector<int> nodes, A_Network *output){return;}
    void isConnected(A_Network network, bool& isConnected){return;}
    void find_neighbours(int node,A_Network network,int distance,vector<int> *neighbours)
    {
      get_unq_neighbors(node, network, distance, neighbours);
      return;
    }

};

#endif
