#ifndef BFS_HPP
#define BFS_HPP

#include "structure_defs.hpp"
#include <iostream>

using namespace std;

// BFS algorithm for A_Network that returns list of visited nodes.
void a_network_bfs(A_Network network, vector<int> &visited) 
{

  // Set up BFS
  visited.clear();
  vector<int> queue;
  queue.push_back(0);
  int head;
  bool bfs_visited;
  int visited_count;

  // Loop through rest of network
  while (!queue.empty()) {
    head = queue.front();
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
	for (int j = 0; j < network.size(); j++) {
	  if (network[head].ListW[i].first == network[j].Row) {
	    queue.push_back(j);
	  }
	}
	//queue.push_back(network[head].ListW[i].first);
      }
      
    }
    visited.push_back(network[head].Row);
  }
  return;
  
}


// DFS algorithm for A_Network that returns list of visited nodes.                                                                                              
void a_network_dfs(A_Network network, vector<int> &visited)
{

  // Set up BFS                                                                                                       
  visited.clear();
  vector<int> queue;
  queue.push_back(0);
  int head;
  bool bfs_visited;
  int visited_count;

  // Loop through rest of network                                                                                 
  while (!queue.empty()) {
    //head = queue.front();
    //queue.erase(queue.begin());
    head = queue.back();
    queue.pop_back();
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
        for (int j = 0; j < network.size(); j++) {
          if (network[head].ListW[i].first == network[j].Row) {
            queue.push_back(j);
          }
        }
        //queue.push_back(network[head].ListW[i].first);                                                      
      }

    }
    visited.push_back(network[head].Row);
  }
  return;

}


// DFS algorithm for A_Network that returns list of visited nodes.                                                                                     
void a_network_dir_dfs(A_Network network, vector<int> &visited)
{

  // Set up BFS                                                                                                                         
  visited.clear();
  vector<int> queue;
  queue.push_back(0);
  int head;

  // Make variables for containing checks
  bool bfs_visited;
  int visited_counter;
  bool is_neighbor;
  int neighbor_counter;

  // Loop through rest of network                                                                           
  while (!queue.empty()) {
    //head = queue.front();                                                                                    
    //queue.erase(queue.begin());                                                                                               
    head = queue.back();
    queue.pop_back();
    int found_neighbor = -1;

    // Search through possible neighbors of the head
    for (int i = 0; i < network.size(); i++) {       
      // Check if network[i] corresponds to an in or out neighbor of network[head]
      is_neighbor = false;
      neighbor_counter = 0;
      while (!is_neighbor && (neighbor_counter < network[i].ListW.size())) {
	if (network[i].ListW[neighbor_counter].first == network[head].Row) {
	  is_neighbor = true;
	  //found_neighbor = i;
	}
	neighbor_counter += 1;
      }
      neighbor_counter = 0;
      while (!is_neighbor && (neighbor_counter < network[head].ListW.size())) {
	if (network[head].ListW[neighbor_counter].first == network[i].Row) {
	  is_neighbor = true;
	  //found_neighbor = i;
	}
	neighbor_counter += 1;
      }
 

      //if (is_neighbor) {
      //	cout << network[i].Row << " is a neighbor of " << network[head].Row << endl;
      //} else {
      //	cout << network[i].Row << " is not a neighbor of " << network[head].Row << endl;
      //}

      if (is_neighbor) {

	// Determine if neighbor is already visited                                                                       
	bfs_visited = false;
	visited_counter = 0;
	while (!bfs_visited && (visited_counter < visited.size())) {
	  if (network[i].Row == visited[visited_counter]) {
	    bfs_visited = true;
	  }
	  visited_counter += 1;
	}
	
	// If neighbor is not already visited, add neighbor to queue
	if (!bfs_visited) {
	  queue.push_back(i);
	}
      }
    }
    //  }
  visited.push_back(network[head].Row);
  }
  return;

}


// Returns the shortest path length to node from every node in network. Path lengths are stored in distance and indexed by the index of network.
void bfs_shortest_paths(int node, A_Network network, vector<int> &distance) 
{

  // Set up distance vector
  distance.clear();
  int node_index;
  distance.resize(network.size());
  for (int i = 0; i < distance.size(); i++) {
    distance[i] = distance.size() + 2;
    if (network[i].Row == node) {
      node_index = i;
    }
  }
  distance[node_index] = 0;

  // Declare needed variables for algorithm
  vector<int> queue;
  queue.push_back(node_index);
  int current_node;
  int neighbor_index;

  // Loop through rest of network
  while (!queue.empty()) {
    current_node = queue.front();
    queue.erase(queue.begin());
    
    // Update system if path to neighbor is longer than the distance to the current node + 1
    for (int neighbor = 0; neighbor < network[current_node].ListW.size(); neighbor++) {

      for (int k = 0; k < network.size(); k++) {
	if (network[k].Row == network[current_node].ListW[neighbor].first) {
	  neighbor_index = k;
	}
      }

      if (distance[neighbor_index] > distance[current_node] + 1) {
	distance[neighbor_index] = distance[current_node] + 1;
	queue.push_back(neighbor_index);
      }

    }
    
  }
  return;

}


// Returns the shortest path length to node from every node in network. Path lengths are stored in distance and indexed by the index of network.        
void dfs_shortest_paths(int node, A_Network network, vector<int> &distance)
{

  // Set up distance vector                                                                   
  distance.clear();
  int node_index;
  distance.resize(network.size());
  for (int i = 0; i < distance.size(); i++) {
    distance[i] = distance.size() + 2;
    if (network[i].Row == node) {
      node_index = i;
    }
  }
  distance[node_index] = 0;

  // Declare needed variables for algorithm                                                                 
  vector<int> queue;
  queue.push_back(node_index);
  int current_node;
  int neighbor_index;

  // Loop through rest of network                                                                    
  while (!queue.empty()) {
    //current_node = queue.front();
    //queue.erase(queue.begin());
    current_node = queue.back();
    queue.pop_back();

    // Update system if path to neighbor is longer than the distance to the current node + 1                                           
    for (int neighbor = 0; neighbor < network[current_node].ListW.size(); neighbor++) {

      for (int k = 0; k < network.size(); k++) {
        if (network[k].Row == network[current_node].ListW[neighbor].first) {
          neighbor_index = k;
        }
      }

      if (distance[neighbor_index] > distance[current_node] + 1) {
        distance[neighbor_index] = distance[current_node] + 1;
        queue.push_back(neighbor_index);
      }

    }

  }
  return;

}


// Returns the shortest path length to node from every node in network. Path lengths are stored in distance and indexed by the index of network.                      
void dir_dfs_shortest_paths(int node, A_Network network, vector<int> &distance)
{

  // Set up distance vector                                                                                                        
  distance.clear();
  int node_index;
  distance.resize(network.size());
  for (int i = 0; i < distance.size(); i++) {
    distance[i] = distance.size() + 2;
    if (network[i].Row == node) {
      node_index = i;
    }
  }
  distance[node_index] = 0;

  // Declare needed variables for algorithm                                                               
  vector<int> queue;
  queue.push_back(node_index);
  int current_node;
  int neighbor_index;

  bool is_neighbor;
  int neighbor_counter;

  // Loop through rest of network                                                   
  while (!queue.empty()) {
    //current_node = queue.front();                                                            
    //queue.erase(queue.begin());                                                            
    current_node = queue.back();
    queue.pop_back();

    // Update system if path to neighbor is longer than the distance to the current node + 1                
    for (int i = 0; i < network.size(); i++) {

      // Check if network[i] corresponds to an in or out neighbor of network[head]                        
      is_neighbor = false;
      neighbor_counter = 0;
      while (!is_neighbor && (neighbor_counter < network[i].ListW.size())) {
        if (network[i].ListW[neighbor_counter].first == network[current_node].Row) {
          is_neighbor = true;
        }
        neighbor_counter += 1;
      }
      neighbor_counter = 0;
      while (!is_neighbor && (neighbor_counter < network[current_node].ListW.size())) {
        if (network[current_node].ListW[neighbor_counter].first == network[i].Row) {
          is_neighbor = true;
        }
        neighbor_counter += 1;
      }

      if ((distance[i] > distance[current_node] + 1) && is_neighbor) {
        distance[i] = distance[current_node] + 1;
        queue.push_back(i);
      }

    }

  }
  return;

}

#endif
