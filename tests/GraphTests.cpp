#include <algorithm>
#include "../include/Graph.hpp"
using namespace ESSENS;

int main()
{
    //Edge petersen[] = {{0,1},{0,4},{0,5},{1,2},{1,6},{2,3},{2,7},{3,4},{3,8},{4,9},{5,7},{5,8},{6,8},{6,9},{7,9}};
    //Graph g(petersen,10,15);
    
    // G17
    Motif motif;
    motif.root = 0;
    motif.nodes = {0,1,2,3,4};
    motif.endnodes = {0,1,2,3,4};
    motif.edges = {{0,1},{0,2},{0,3},{1,2},{1,3},{1,4}};
    
    printMotif(motif);

    std::cout << "Motif Adj Matrix" << std::endl;

    // Edge e;

    for (Node u: motif.nodes)
    {
        for (Node v: motif.nodes)
        {
            int found = 0;
            for (Edge e: motif.edges)
            {
                if (e.u == u && e.v == v)
                {
                    found = 1;
                    break;
                } 
            }
            std::cout << found << " ";
        }
        std::cout << std::endl;
    }
    
    std::map<std::pair<Node,Node>,int> distances;
    for (Node u: motif.nodes)
    {
        for (Node v: motif.nodes)
        {
            distances.insert({{u,v},-1});
        }
    }
    all_pair_shortest_lengths(motif, distances);

    std::cout << "Final counts" << std::endl;
    // for (auto node_pair: distances)
    //     std::cout << node_pair.first.first << " " << node_pair.first.second << ": " << node_pair.second << std::endl; 

    for (Node u: motif.nodes)
    {
        for (Node v: motif.nodes)
        {
            std::cout << distances[{u,v}] << " ";
        }
        std::cout << std::endl;
    }
}
