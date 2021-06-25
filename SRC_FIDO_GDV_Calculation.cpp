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

Get_all_combinations_with_node_at_left_end(A_Network graph,int node,vector<Combination> *list_of_combinations_with_node)
{

}
void GDV_vector_calculation(A_Network graph,vector<GDVMetric>* graph_GDV,  vector<OrbitMetric> orbits,int p)
{
    vector<Combination> list_of_combinations;
    for(int i=0; i<graph.size(); i++){
        
        vector<Combination> list_of_combinations_with_node;
        ADJ_Bundle node = graph[i];
        Get_all_combinations_with_node_at_left_end(graph,node.Row,list_of_combinations_with_node);
    }
    Get_all_back_edges_with_combinations(list_of_combinations_with_node);
}