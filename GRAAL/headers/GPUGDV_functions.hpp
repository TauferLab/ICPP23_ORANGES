// Required directives are to be added from ESSENS to use the functionality 

#ifndef GPUGDV_FUNCTIONS_HPP
#define GPUGDV_FUNCTIONS_HPP

#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
#include "ADJ/network_defs.hpp"
#include "ADJ/create_network.hpp"
#include "combinations_raw.hpp"
#include "class_definitions.hpp"
//#include "ADJ/traversal.hpp"
//#include <iostream>
#include "stdio.h"
//#include "bfs.hpp"
#include "print_disconnected_graph.hpp"
#include "raw_vecs.hpp"
#include "gpubfs.hpp"


//#include "/home/pnbell/Src_GraphAlignment/ESSENS-master/Core/Basic_IO/Format/Level0/structure_defs.hpp"                                                                                         
using namespace std;

#pragma once

class GPUGDV_functions{

public:
       
  int fact(int n){
    return (n==0) || (n==1) ? 1 : n*fact(n-1);
  }

  void find_combinations(int set[], int n, int size, intvecvec& output, umpire::Allocator &alloc){
     Combinations_raw c;
     c.getCombination(set, n, size, output);
  }

  // Vector nodes should correspond to actual node labels, not node indices.
  FIDO_HOST_DEVICE void inducedSubgraph_raw(A_Network_raw &network, intvec nodes, A_Network_raw& output)
  {

    // Set up output for new subgraph. 

    // The induced subgraph should probably be cleared outside of here,
    // and its sublists freed up or recleared. For now, though, just clearing...
    clear_adjlist(output);  

    // Create needed objects for function.                                                                                                                                                                                                   
    int subgraph_nodes = nodes.veclen;
    edgevec neighbor_list;
    Adjlist temp_bundle;
    Edge temp_edge;

    // Check the neighbors of the subgraph nodes to determine where edges need to be.                                                                                                                                                        
    for (int i = 0; i < subgraph_nodes; i++) {

      // Find which node in the network corresponds to nodes[i] and store it's neighbors.
      bool node_found = false;
      int graph_count = 0;
      while (!node_found) {
	if (network.vec[graph_count].Row == nodes.vec[i]) {
	  node_found = true;
	  neighbor_list.vec = network.vec[graph_count].ListW;
	}
	graph_count += 1;
      }

      // Ensure subgraph contains next node.              
      //temp_bundle.Row = nodes.vec[i];
      pushback_adjlist(output, new_adjlist(subgraph_nodes, 0));
      output.vec[output.veclen-1].Row=nodes.vec[i];
      

      // Check which of the nodes[j] correspond to neighbors of nodes[i].
      for (int j = 0; j < subgraph_nodes; j++) {
	for (int k = 0; k < neighbor_list.veclen; k++) {
	  if (neighbor_list.vec[k].first == nodes.vec[j]) {
	    temp_edge.first = nodes.vec[j];
	    temp_edge.second = neighbor_list.vec[k].second;
            pushback_edgevec(output.vec[i].ListW, temp_edge);
	    //output[i].ListW.push_back(temp_edge);
	  }
	}
      }

    }

    return;
  };

  // Evaluate whether network is a connected graph.  Store result in isConnected.                  
  FIDO_HOST_DEVICE void isConnected_raw(A_Network_raw &network, bool& isConnected)
  {
    intvec visited_nodes = new_intvec(network.nodes_len);

    if (network.nodes_len != 0) {
      a_network_dir_dfs_raw(network, visited_nodes);
      
      // Make connectedness check
      if (visited_nodes.veclen == network.nodes_len) {
	isConnected = true;
      }
      else {
	isConnected = false;
      }
    } else {
      printf("Error in function isConnected: input variable network is empty.\n");
    }
    delete_intvec(visited_nodes);
    return;
  };


  // Calculate degree signature for network
  FIDO_HOST_DEVICE void degree_signature_raw(A_Network_raw &network, intvec &deg_sig)
  {
    if (network.nodes_len != 0) {
      clear_intvec(deg_sig);
      for (int i = 0; i < network.nodes_len; i++) {
        pushback_intvec(deg_sig, network[i].ListW.veclen);
      }
    } else {
      printf("Error in function degree_signature: input variable network is empty.\n");
    }

    return;
  }

  // Calculate distance signature for network
  FIDO_HOST_DEVICE void distance_signature_raw(int node, A_Network_raw &network, intvec &dist_sig)
  {
    if (network.nodes_len != 0) {

      // Start by calculating shortest paths in graph.
      
      // This could be externally allocated
      intvec shortest_paths = new_intvec(network.nodes_len);      

      //vector<int> shortest_paths;
      dir_dfs_shortest_paths_raw(node, network, shortest_paths);
	
      // Then use the shortest paths to calculate the distance signature.
      clear_intvec(dist_sig);
      
      for (int i = 0; i < shortest_paths.veclen; i++) {
	dist_sig.vec[shortest_paths.vec[i]] += 1;
      }

      delete_intvec(shortest_paths);      

    } else {
      printf("Error in function distance_signature: input variable network is empty.\n");
    }
    return;
  }

  FIDO_HOST_DEVICE void orbit_filter_2(orbvec& orbits, orbvec& filtered_orbits) {
      //filtered_orbits.push_back(orbits[0]);
      pushback_orbvec(filtered_orbits, orbits.vec[0]);
  }

  FIDO_HOST_DEVICE void orbit_filter_3(orbvec& orbits, orbvec& filtered_orbits) {
      //filtered_orbits.push_back(orbits[1]);
      pushback_orbvec(filtered_orbits, orbits.vec[1]);
      //filtered_orbits.push_back(orbits[2]);
      pushback_orbvec(filtered_orbits, orbits.vec[2]);
      //filtered_orbits.push_back(orbits[3]);
      pushback_orbvec(filtered_orbits, orbits.vec[3]);
  }

  FIDO_HOST_DEVICE void orbit_filter_4(orbvec& orbits, orbvec& filtered_orbits) {
      //filtered_orbits.push_back(orbits[4]);
      pushback_orbvec(filtered_orbits, orbits.vec[4]);
      //filtered_orbits.push_back(orbits[5]);
      pushback_orbvec(filtered_orbits, orbits.vec[5]);
      //filtered_orbits.push_back(orbits[6]);
      pushback_orbvec(filtered_orbits, orbits.vec[6]);
      //filtered_orbits.push_back(orbits[7]);
      pushback_orbvec(filtered_orbits, orbits.vec[7]);
      //filtered_orbits.push_back(orbits[8]);
      pushback_orbvec(filtered_orbits, orbits.vec[8]);
      //filtered_orbits.push_back(orbits[9]);
      pushback_orbvec(filtered_orbits, orbits.vec[9]);
      //filtered_orbits.push_back(orbits[10]);
      pushback_orbvec(filtered_orbits, orbits.vec[10]);
      //filtered_orbits.push_back(orbits[11]);
      pushback_orbvec(filtered_orbits, orbits.vec[11]);
      //filtered_orbits.push_back(orbits[12]);
      pushback_orbvec(filtered_orbits, orbits.vec[12]);
      //filtered_orbits.push_back(orbits[13]);
      pushback_orbvec(filtered_orbits, orbits.vec[13]);
      //filtered_orbits.push_back(orbits[14]);
      pushback_orbvec(filtered_orbits, orbits.vec[14]);
  }

  FIDO_HOST_DEVICE void orbit_filter_5(orbvec& orbits, orbvec& filtered_orbits) {
      //filtered_orbits.push_back(orbits[15]);
      pushback_orbvec(filtered_orbits, orbits.vec[15]);
      //filtered_orbits.push_back(orbits[16]);
      pushback_orbvec(filtered_orbits, orbits.vec[16]);
      //filtered_orbits.push_back(orbits[17]);
      pushback_orbvec(filtered_orbits, orbits.vec[17]);
      //filtered_orbits.push_back(orbits[18]);
      pushback_orbvec(filtered_orbits, orbits.vec[18]);
      //filtered_orbits.push_back(orbits[19]);
      pushback_orbvec(filtered_orbits, orbits.vec[19]);
      //filtered_orbits.push_back(orbits[20]);
      pushback_orbvec(filtered_orbits, orbits.vec[20]);
      //filtered_orbits.push_back(orbits[21]);
      pushback_orbvec(filtered_orbits, orbits.vec[21]);
  }

  FIDO_HOST_DEVICE void orbit_filter_default(orbvec& orbits, orbvec& filtered_orbits) {
      //filtered_orbits.push_back(orbits[0]);
      pushback_orbvec(filtered_orbits, orbits.vec[0]);
      //filtered_orbits.push_back(orbits[1]);
      pushback_orbvec(filtered_orbits, orbits.vec[1]);
      //filtered_orbits.push_back(orbits[2]);
      pushback_orbvec(filtered_orbits, orbits.vec[2]);
      //filtered_orbits.push_back(orbits[3]);
      pushback_orbvec(filtered_orbits, orbits.vec[3]);
      //filtered_orbits.push_back(orbits[4]);
      pushback_orbvec(filtered_orbits, orbits.vec[4]);
      //filtered_orbits.push_back(orbits[5]);
      pushback_orbvec(filtered_orbits, orbits.vec[5]);
      //filtered_orbits.push_back(orbits[6]);
      pushback_orbvec(filtered_orbits, orbits.vec[6]);
      //filtered_orbits.push_back(orbits[7]);
      pushback_orbvec(filtered_orbits, orbits.vec[7]);
      //filtered_orbits.push_back(orbits[8]);
      pushback_orbvec(filtered_orbits, orbits.vec[8]);
      //filtered_orbits.push_back(orbits[9]);
      pushback_orbvec(filtered_orbits, orbits.vec[9]);
      //filtered_orbits.push_back(orbits[10]);
      pushback_orbvec(filtered_orbits, orbits.vec[10]);
      //filtered_orbits.push_back(orbits[11]);
      pushback_orbvec(filtered_orbits, orbits.vec[11]);
      //filtered_orbits.push_back(orbits[12]);
      pushback_orbvec(filtered_orbits, orbits.vec[12]);
      //filtered_orbits.push_back(orbits[13]);
      pushback_orbvec(filtered_orbits, orbits.vec[13]);
      //filtered_orbits.push_back(orbits[14]);
      pushback_orbvec(filtered_orbits, orbits.vec[14]);
      //filtered_orbits.push_back(orbits[15]);
      pushback_orbvec(filtered_orbits, orbits.vec[15]);
      //filtered_orbits.push_back(orbits[16]);
      pushback_orbvec(filtered_orbits, orbits.vec[16]);
      //filtered_orbits.push_back(orbits[17]);
      pushback_orbvec(filtered_orbits, orbits.vec[17]);
      //filtered_orbits.push_back(orbits[18]);
      pushback_orbvec(filtered_orbits, orbits.vec[18]);
      //filtered_orbits.push_back(orbits[19]);
      pushback_orbvec(filtered_orbits, orbits.vec[19]);
      //filtered_orbits.push_back(orbits[20]);
      pushback_orbvec(filtered_orbits, orbits.vec[20]);
      //filtered_orbits.push_back(orbits[21]);
      pushback_orbvec(filtered_orbits, orbits.vec[21]);
  }

  FIDO_CONSTANT void (*subfunc_ptrs[5]) = {
      orbit_filter_default,
      orbit_filter_2,
      orbit_filter_3,
      orbit_filter_4,
      orbit_filter_5
  };

  FIDO_HOST_DEVICE void orbit_filter_raw(orbvec& orbits, int nodes, orbvec& filtered_orbits)
  {
      if (nodes < 2 || nodes > 5) {
          nodes = 1;
      }
      (*subfunc_ptrs[nodes-1])(orbits, filtered_orbits);
  }

};


#endif


