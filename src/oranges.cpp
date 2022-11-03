/**
 * @file oranges.cpp
 * @author Ali Khan (alikhan@my.unt.edu)
 * @brief 
 * @version 0.1
 * @date 2022-09-13
 * 
 * @copyright Copyright (c) 2022
 * 
 */
// #include "structure_defs.hpp"
// #include "input_to_network.hpp"
// #include "printout_others.hpp"
// #include "printout_network.hpp"


// // #include "ADJ/find_Xneighbors.hpp"


// //INPUT HEADERS
// #include "translate_from_input.hpp"
// #include "input_to_network.hpp"
// #include"structure_defs.hpp"

// //OUTPUT HEADERS
// #include "printout_network.hpp"
// #include "printout_others.hpp"

// // // HEADERS
// // #include "CMPB_structure.hpp"
// // #include  "subgraph_extract.hpp"
// // #include "createN_Struct.hpp"

// #include "ADJ/degree_centrality.hpp"

// #include "add_multiple_edge.hpp"
// #include "ADJ/create_network.hpp"
// /*** All Headers Required From ESSENS **/

// #include "reorder.hpp"

#include <time.h>
#include <math.h>
#include <omp.h>
#include <algorithm>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <utility>

const unsigned int MAX_NODES = 5;
const unsigned int MAX_ORBIT = 73;

typedef unsigned int Node;
// class Node
// {
//     int value;
//     bool isEndnode = false;
// };

class Edge
{
    public:
    Node u;
    Node v;
    Edge(Node _u, Node _v)
    {
        u = _u;
        v = _v; 
    }
};

class Graph
{
    public:
    int N;
    int M;
    std::vector<Node> Nodes;
    std::map<Node,std::vector<Node>> adjList;
    void addEdge(Edge e);
    std::vector<Node> neighbors(Node u);
};

class Motif
{
    public:
    Node root; // Node is basically int. Edge is basically pair<Node,Node>
    //std::vector<std::pair<Node,bool>> nodes; // pair with first int the node id and the second whether the node is a leaf
    std::map<Node,bool> nodes;
    std::vector<Edge> edges;
};


void combinations(std::vector<Edge> &edges, int num_neighbors, std::vector<std::vector<Edge>> &combos)
{
    // enumerate combinations of edges
    return;
}

void combinations(std::vector<Node> &nodes, int num_neighbors, std::vector<std::vector<Node>> &combos)
{
    // enumerate combinations of nodes
    return;
}

int getOrbit(Motif motif, Node root)
{
    return;
}

/**
 * @brief Compute the graphlet degree vector (GDV) for each node
 * 
 * @param g a 
 */
void computeGdv(Graph &g, std::map<Node,unsigned int[MAX_ORBIT]> gdv) // std::map<Node,std::vector<unsigned int>> gdv)
{
    std::queue<Motif> subtreeQueue;
    Node root;
    std::vector<std::vector<Node>> leaf_combinations;
    std::vector<std::vector<Edge>> edge_combinations;
    Motif tree;
    Motif newTree;
    std::set<Edge> edge_set;
    std::vector<Node> neighbors;
    for (std::vector<Node>::iterator it = g.Nodes.begin(); it != g.Nodes.end(); ++it)
    {
        root = *it;

        // get neighbors
        std::vector<Node> neighbors;
        for (Node u: g.neighbors(root))
        {
            if (u > root)
            {
                neighbors.push_back(u);
            }
        }
        
        //std::vector<std::vector<Node>> leaf_combinations;
        combinations(neighbors, MAX_NODES, leaf_combinations);

        // initialize the subtree queue
        for (std::vector<Node> combo : leaf_combinations)
        {
            // create motif
            Motif tree;
            tree.root = root;
            for (Node w: combo)
            {
                if (w > root)
                {
                    tree.nodes.insert(std::make_pair(w,true)); // endnode = True
                    tree.edges.push_back(Edge(w,root));
                }
            }

            subtreeQueue.push(tree);

            //orbitCounts[getOrbit(tree,root)] += 1;
            gdv[root][getOrbit(tree,root)] += 1;
        }

        // main loop
        while (!subtreeQueue.empty())
        {
            // Motif tree = subtreeQueue.front();
            subtreeQueue.pop();
            if (tree.nodes.size() <= MAX_NODES)
            {
                std::vector<Edge> neighboringEdges;
                for (auto const& nodepair: tree.nodes)
                {
                    Node v = nodepair.first;
                    bool v_isEndnode = nodepair.second;
                    for (Node w : g.neighbors(v))
                    {
                        if (v_isEndnode && w > root && tree.nodes.find(w) == tree.nodes.end())
                        {
                            neighboringEdges.push_back(Edge(w,v));
                        }
                    }
                }

                //std::vector<std::vector<Edge>> edge_combinations;
                combinations(neighboringEdges, MAX_NODES - tree.nodes.size() + 1, edge_combinations);

                for (std::vector<Edge> combo : edge_combinations)
                {
                    // filter out combinations that include multiple of the same edge
                    // std::set<Edge> edge_set;
                    bool skip = false;
                    for (Edge e : combo)
                    {
                        if (edge_set.find(e) != edge_set.end())
                        {
                            edge_set.insert(e);
                        }
                        else 
                        {
                            skip = true;
                            break;
                        }
                    }
                    if (skip)
                    {
                        continue;
                    }

                    //create new subtree
                    // Motif newTree;
                    newTree.root = tree.root;
                    // copy edges
                    for (int i=0; i<tree.edges.size(); i++)
                    {
                        newTree.edges.push_back(tree.edges[i]); 
                    }
                    // copy nodes
                    for (auto const& elem: tree.nodes)
                    {
                        newTree.nodes.insert(elem); 
                    }
                    // add new edges
                    for (Edge e : combo)
                    {
                        newTree.edges.push_back(e);
                        // unmark endnodes
                        newTree.nodes[e.v] = false;
                        
                        // add new node (the second node in the edge) and mark it as an endnode
                        newTree.nodes.insert(std::make_pair(e.u,true));
                        
                        // remove endnode designation
                        newTree.nodes[e.u] = false;
                    }

                    subtreeQueue.push(newTree);

                    gdv[root][getOrbit(tree,root)] += 1;
                }
            }

        }
    }
}
