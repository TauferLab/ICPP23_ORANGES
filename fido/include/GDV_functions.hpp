// Required directives are to be added from ESSENS to use the functionality 

#ifndef GDV_FUNCTIONS_HPP
#define GDV_FUNCTIONS_HPP

//#include "/home/pnbell/Src_GraphAlignment/ESSENS-master/Core/Basic_IO/Format/Level0/structure_defs.hpp"                                                                                         
#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
#include "ADJ/network_defs.hpp"
#include "ADJ/create_network.hpp"
#include "combinations.hpp"
#include "class_definitions.hpp"
//#include "ADJ/traversal.hpp"
#include <iostream>
#include "bfs.hpp"
#include "print_disconnected_graph.hpp"

using namespace std;

#pragma once

class GDV_functions{

public:

  // List all combinations of size from set.  Store in output.                     
  void find_combinations ( int set[], int n,int size, vector<vector<int>> *output){
    Combinations c;
    c.getCombination(set,n,size,output);
    return;
  }

  // Create a graph output using nodes from network.                                                                                                                             
  // Vector nodes should correspond to actual node labels, not node indices.
  void inducedSubgraph(A_Network network, const vector<int> &nodes, A_Network& output)
  {

    // Set up output for new subgraph.                                                                                                                                                                                                   
    output.clear();

    // Create needed objects for function.                                                                                                                                                                                                   
    int subgraph_nodes = nodes.size();
    vector<int_double> neighbor_list;
    ADJ_Bundle temp_bundle;
    int_double temp_edge;

    // Check the neighbors of the subgraph nodes to determine where edges need to be.                                                                                                                                                        
    for (int i = 0; i < subgraph_nodes; i++) {

      // Find which node in the network corresponds to nodes[i] and store it's neighbors.
      bool node_found = false;
      int graph_count = 0;
      while (!node_found) {
	      if (network[graph_count].Row == nodes[i]) {
	        node_found = true;
	        neighbor_list = network[graph_count].ListW;
	        //cout << nodes[i] << ": ";
	        //for (int j = 0; j < neighbor_list.size(); j++) {
	        //  cout << neighbor_list[j].first << ", ";
	        //}
	        //cout << endl;
	      }
	      graph_count += 1;
      }

      // Ensure subgraph contains next node.              
      temp_bundle.Row = nodes[i];
      output.push_back(temp_bundle);

      // Check which of the nodes[j] correspond to neighbors of nodes[i].
      for (int j = 0; j < subgraph_nodes; j++) {
	      for (int k = 0; k < neighbor_list.size(); k++) {
	        if (neighbor_list[k].first == nodes[j]) {
	          temp_edge.first = nodes[j];
	          temp_edge.second = neighbor_list[k].second;
	          output[i].ListW.push_back(temp_edge);
	        }
	      }
      }

    }

    return;
  };


  // Evaluate whether network is a connected graph.  Store result in isConnected.                  
  void isConnected(A_Network network, bool& isConnected)
  {
    A_Network spanning_tree;
    vector<int> visited_nodes;

    if (!network.empty()) {
      a_network_dir_dfs(network, visited_nodes);
      
      // Make connectedness check
      if (visited_nodes.size() == network.size()) {
	isConnected = true;
      }
      else {
	isConnected = false;
      }
    } else {
      cout << "Error in function isConnected: input variable network is empty." << endl;
    }

    return;
  };

  // List neighbors up to distance from node in network.                                                               
  void find_neighbours(int node,A_Network network,int distance,vector<int> *neighbours)
  {
    get_unq_neighbors(node, network, distance, neighbours);
    sort(neighbours->begin(), neighbours->begin()+neighbours->size());
    unsigned int index = 0;
    for(index; index<neighbours->size(); index++) {
      if((*neighbours)[index] > node)
        break;
    }
    neighbours->erase(neighbours->begin(), neighbours->begin()+index);
    return;
  }

  // Calculate degree signature for network
  void degree_signature(A_Network network, vector<int> &deg_sig)
  {
    if (!network.empty()) {
      deg_sig.clear();
      for (int i = 0; i < network.size(); i++) {
	deg_sig.push_back(network[i].ListW.size());
      }
    } else {
      cout << "Error in function degree_signature: input variable network is empty." << endl;
    }

    return;
  }

  // Calculate distance signature for network
  void distance_signature(int node, A_Network network, vector<int> &dist_sig)
  {
    if (!network.empty()) {

      // Start by calculating shortest paths in graph.
      vector<int> shortest_paths;
      dir_dfs_shortest_paths(node, network, shortest_paths);
	
      // Then use the shortest paths to calculate the distance signature.
      dist_sig.clear();
      dist_sig.resize(6, 0);
      for (int i = 0; i < shortest_paths.size(); i++) {
	      dist_sig[shortest_paths[i]] += 1;
      }

    } else {
      cout << "Error in function distance_signature: input variable network is empty." << endl;
    }
    return;
  }

  vector<OrbitMetric> orbit_filter( vector<OrbitMetric>& orbits, int nodes, vector<OrbitMetric>& filtered_orbits)
  {
     switch(nodes) {
      case 2 :
          filtered_orbits.push_back(orbits[0]);
          break;
      case 3 :
          filtered_orbits.push_back(orbits[1]);
          filtered_orbits.push_back(orbits[2]);
          filtered_orbits.push_back(orbits[3]);
         break;
      case 4 :
            filtered_orbits.push_back(orbits[4]);
            filtered_orbits.push_back(orbits[5]);
            filtered_orbits.push_back(orbits[6]);
            filtered_orbits.push_back(orbits[7]);
            filtered_orbits.push_back(orbits[8]);
            filtered_orbits.push_back(orbits[9]);
            filtered_orbits.push_back(orbits[10]);
            filtered_orbits.push_back(orbits[11]);
            filtered_orbits.push_back(orbits[12]);
            filtered_orbits.push_back(orbits[13]);
            filtered_orbits.push_back(orbits[14]);
         break;
      case 5 :
            filtered_orbits.push_back(orbits[15]);
            filtered_orbits.push_back(orbits[16]);
            filtered_orbits.push_back(orbits[17]);
            filtered_orbits.push_back(orbits[18]);
            filtered_orbits.push_back(orbits[19]);
            filtered_orbits.push_back(orbits[20]);
            filtered_orbits.push_back(orbits[21]);
         break;
      default :
            filtered_orbits.push_back(orbits[0]);
            filtered_orbits.push_back(orbits[1]);
            filtered_orbits.push_back(orbits[2]);
            filtered_orbits.push_back(orbits[3]);
            filtered_orbits.push_back(orbits[4]);
            filtered_orbits.push_back(orbits[5]);
            filtered_orbits.push_back(orbits[6]);
            filtered_orbits.push_back(orbits[7]);
            filtered_orbits.push_back(orbits[8]);
            filtered_orbits.push_back(orbits[9]);
            filtered_orbits.push_back(orbits[10]);
            filtered_orbits.push_back(orbits[11]);
            filtered_orbits.push_back(orbits[12]);
            filtered_orbits.push_back(orbits[13]);
            filtered_orbits.push_back(orbits[14]);
            filtered_orbits.push_back(orbits[15]);
            filtered_orbits.push_back(orbits[16]);
            filtered_orbits.push_back(orbits[17]);
            filtered_orbits.push_back(orbits[18]);
            filtered_orbits.push_back(orbits[19]);
            filtered_orbits.push_back(orbits[20]);
            filtered_orbits.push_back(orbits[21]);
          break;
   }
    return filtered_orbits;
  }

};


#endif


