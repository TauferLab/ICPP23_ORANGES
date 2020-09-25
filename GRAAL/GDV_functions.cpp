// Required directives are to be added from ESSENS to use the functionality 

#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
#include "GDV_functions.hpp"
using namespace std;

GDV_functions::find_combination ( int set[], int size, vector<vector<int>>& output)
{
    return;
};
GDV_functions::inducedSubGraph(A_Network network, vector<int> nodes, A_Network& output)
{
    return;
}
GDV_functions::isConnected(A_Network network, bool& isConnected)
{
    return;
}

GDV_functions::find_neighbours(int node,A_Network network,int distance,vector<int> &neighbours)
{
    get_unq_neighbors(node,network,k,&neighbours);
    return;
}

