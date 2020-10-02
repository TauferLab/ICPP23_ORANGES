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
  queue.push_back(network[0].Row);
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
	queue.push_back(network[head].ListW[i].first);
      }

    }
    visited.push_back(head);
  }
  return;

}


void bfs_shortest_paths(int node, A_Network network, vector<int> &distance) 
{

  // Set up distance vector
  distance.clear();
  distance.resize(network.size());
  for (int i = 0; i < distance.size(); i++) {
    distance[i] = distance.size() + 2;
  }
  distance[node] = 0;

  // Declare needed variables for algorithm
  vector<int> queue;
  queue.push_back(node);
  int current_node;

  // Loop through rest of network
  while (!queue.empty()) {
    current_node = queue.front();
    queue.erase(queue.begin());
    
    // Update system if path to neighbor is longer than the distance to the current node + 1
    for (int neighbor = 0; neighbor < network[current_node].ListW.size(); neighbor++) {
      if (distance[network[current_node].ListW[neighbor].first] > distance[current_node] + 1) {
	distance[network[current_node].ListW[neighbor].first] = distance[current_node] + 1;
	queue.push_back(network[current_node].ListW[neighbor].first);
      }
    }
    
  }
  return;

}



#endif
