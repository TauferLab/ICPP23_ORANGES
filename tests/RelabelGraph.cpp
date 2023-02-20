// from chatgpt

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include "../include/graph.hpp"

using namespace ESSENS;

#if _WIN32  //win32 or win64
    #define INPUT_GRAPH_PATH ".\\data\\karate_undir.txt"
#else  //linux
    #define INPUT_GRAPH_PATH "./data/karate_undir.txt"
#endif

int main() {
    // Edge karate[] = {
    //     {0,9},{0,14},{0,15},{0,16},{0,19},{0,20},{0,21},{0,23},{0,24},{0,27},{0,28},{0,29},{0,30},{0,31},{0,32},{0,33},{2,1},{3,1},{3,2},{4,1},{4,2},{4,3},{5,1},{6,1},{7,1},{7,5},{7,6},{8,1},{8,2},{8,3},{8,4},{9,1},{9,3},{10,3},{11,1},{11,5},{11,6},{12,1},{13,1},{13,4},{14,1},{14,2},{14,3},{14,4},{17,6},{17,7},{18,1},{18,2},{20,1},{20,2},{22,1},{22,2},{26,24},{26,25},{28,3},{28,24},{28,25},{29,3},{30,24},{30,27},{31,2},{31,9},{32,1},{32,25},{32,26},{32,29},{33,3},{33,9},{33,15},{33,16},{33,19},{33,21},{33,23},{33,24},{33,30},{33,31},{33,32}
    // };
    // Graph graph(karate,34,77);
    std::vector<Edge> edgelist = read_graph(INPUT_GRAPH_PATH);
    Graph graph(edgelist);

    std::vector<Node> new_ordering;
    kcore_ordering(graph,new_ordering);
    std::cout << "old:new" << std::endl;
    for (size_t i = 0; i < new_ordering.size(); i++)
    {
        std::cout << i << ":" << new_ordering[i] << std::endl;
    }
    graph.relabel_graph(new_ordering);
    graph.printGraph();
    return 0;
}
