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
//#include "ADJ/traversal.hpp"
#include <iostream>

using namespace std;

#pragma once

class GDV_functions{

public:

  // List all combinations of size from set.  Store in output.                                                                                                                                                                                
  void find_combinations(int set[], int size, vector<vector<int>>& output)
  {
    
  };

  // Create a graph output using nodes from network.                                                                                                                                                                                          
  void inducedSubgraph(A_Network network, vector<int> nodes, A_Network& output)
  {

    // Store number of nodes being passed.                                                                                                                                                                                                   
    int subgraph_nodes = nodes.size();
    A_Network* subgraph = &output;
    output.clear();

    // Create needed objects for function.                                                                                                                                                                                                   
    vector<int> nodei;
    nodei.push_back(0);
    vector<int> neighbor_list;

    // Variables for storing ADJ information                  
    ADJ_Bundle temp_bundle;
    int_double temp_edge;

    // Check the neighbors of the subgraph nodes to determine where edges need to be.                                                                                                                                                        
    for (int i = 0; i < subgraph_nodes; i++) {

      nodei[0] = network[nodes[i]].Row;
      get_neighbors(nodei, network, &neighbor_list);

      // Ensure subgraph contains next node.              
      temp_bundle.Row = i;
      output.push_back(temp_bundle);

      for (int j = 0; j < subgraph_nodes; j++) {

        // Retrieve edge weight for potential edge                
        for (int k = 0; k < network[nodes[i]].ListW.size(); k++) {
          if (network[nodes[i]].ListW[k].first == network[nodes[j]].Row) {
            temp_edge.first = j;
            temp_edge.second = network[nodes[i]].ListW[k].second;
          }
        }

        // Determine if j is a neighbor of i by looping through the neighbors of i                
        for (int k = 0; k < neighbor_list.size(); k++) {
          if (neighbor_list[k] == network[nodes[j]].Row) {
            output[i].ListW.push_back(temp_edge);
          }
        }

      }

    }

    return;
  };

  // BFS algorithm for A_Network that returns list of visited nodes.
  void a_network_bfs(A_Network network, vector<int> &visited) {

    // Set up BFS
    visited.clear();
    vector<int> queue;
    queue.push_back(network[0].Row);
    int head;
    bool bfs_visited;
    int visited_count;

    // Loop through rest of network
    while (!queue.empty()) {
      head = queue[0];
      queue.erase(queue.begin());
      for (int i = 0; i < network[head].ListW.size(); i++) {
	
	// Determine if neighbor is already visited
	bfs_visited = false;
	visited_count = 0;
	while (!bfs_visited && visited_count < visited.size()) {
	  if (network[head].ListW[i].first == visited[visited_count]) {
	    bfs_visited = true;
	  } 
	  visited_count += 1;
	}

	// If neighbor not already visted, add neighbor to queue
	if (!bfs_visited) {
	  queue.push_back(network[head].ListW[i].first);
	}

      }
      visited.push_back(head);
    }
    return;

  };

  // Evaluate whether network is a connected graph.  Store result in isConnected.                  
  void isConnected(A_Network network, bool& isConnected)
  {
    A_Network spanning_tree;
    vector<int> visited_nodes;
    a_network_bfs(network, visited_nodes);
    
    // Make connectedness check
    if (visited_nodes.size() == network.size()) {
      isConnected = true;
    }
    else {
      isConnected = false;
    }

    return;
  };

  // List neighbors up to distance from node in network.                                                               
  void find_neighbours(int node,A_Network network,int distance,vector<int> &neighbours)
  {
    get_unq_neighbors(node, network, distance, &neighbours);
    return;
  };

  // Can be used for printing disconnected graphs.             
  void print_disconnected_network(A_Network &network) 
  {
    for (int i = 0; i < network.size(); i++) {
      cout << "-------------------------------------" << endl;
      cout << network[i].Row << ":" << endl;
      for (int j = 0; j < network[i].ListW.size(); j++) {
	cout << network[i].ListW[j].first << " , " << endl;
      }
      cout << "-------------------------------------" << endl;
    }
    return;
  };

};


#endif
