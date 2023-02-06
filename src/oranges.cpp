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


#include <time.h>
#include <math.h>
#include <omp.h>
#include <algorithm>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <utility>

#include "../include/graph.hpp"
#include "../include/orbitmatch.hpp"
#include "../include/combinations.hpp"

const unsigned int MAX_NODES = 5;
const unsigned int MAX_ORBIT = 72;

using namespace ESSENS;


// for each node ine motif, it takes the orbit found by the orbitmatcher and updates the gdv accordingly
bool updateGDV(const Motif &motif, const OrbitMatcher orbitmatcher, std::vector<std::array<int,MAX_ORBIT>> &gdv)
{
    int orbit;
    size_t i;
    Node node;
    bool success = true;
    std::vector<int> orbits = orbitmatcher.get_orbits(motif);

    for (i = 0; i < motif.nodes.size(); i++)
    {
        orbit = orbits[i];
        node = motif.nodes[i];
        if (orbit != -1)
        {
            gdv[node][orbit]++;
        } else{
            success = false;
        }
    }
    return success;
}

// creates new motifs using backedges of the tree and updates the gdv accordingly
void update_gdv_backedges(Motif const &tree, Graph const &g, std::vector<std::array<int,MAX_ORBIT>> &gdv)
{
    // 1. get backedges
    // 2. create new motifs
    // 3. update gdv
}

// skip this combo if it has any backedges ie {(1,2), (1,3)}
// prevents repeating work on the same tree
bool skip_combo(Motif &motif, std::vector<int> &edge_combination, std::vector<Edge> &neighboringEdges)
{
    Node node;
    std::unordered_set<Node> node_set;
    for (auto edge_idx : edge_combination)
    {
        node = neighboringEdges[edge_idx].u;
        if (node_set.find(node) == node_set.end())
        {
            node_set.insert(node);
        } else
        {
            // skip this combo since it is a backedge; another edge combination will cover it
            return true;
        }
    }
    return false;
}

/**
 * @brief Compute the graphlet degree vector (GDV) for each node
 * 
 * @param g a 
 */
void computeGDV(Graph &g, const OrbitMatcher &orbitmatcher, std::vector<std::array<int,MAX_ORBIT>> &gdv) // std::map<Node,std::vector<unsigned int>> gdv)
{
    std::queue<Motif> subtreeQueue;
    Edge e;
    
    std::vector<Node> neighbors;
    std::vector<Edge> neighboringEdges;

    for (Node root: g.Nodes)
    {
        std::cout << "\nRoot=" << root << std::endl;

        // get neighbors
        neighbors.clear();
        for (Node u: g.neighbors(root))
        {
            // we only want neighbors that have higher value than the root
            if (u > root)
            {
                neighbors.push_back(u);
            }
        }
        std::cout << "neighbors: [";
        for (auto neighbor: neighbors)
            std::cout << neighbor << ",";
        std::cout << "]" << std::endl;

        // if no valid neighbors, skip to new root
        if (neighbors.size() == 0)
        {
            continue;
        }

        // setup combination generator. n choose 1..k
        AllCombinationGenerator combo_gen(neighbors.size(), std::min((int) neighbors.size(),(int)MAX_NODES - 1));
        
        // initialize the subtree queue
        while(!combo_gen.done)
        {
            std::cout << "combo " << combo_gen.combo_cnt << ": [";
            for (auto idx: combo_gen.indices)
                std::cout << neighbors[idx] << ",";
            std::cout << "]" << std::endl;

            // create motif
            Motif tree;
            tree.root = root;
            tree.nodes.push_back(root);
            for (int idx: combo_gen.indices)
            {
                tree.nodes.push_back(neighbors[idx]);
                tree.endnodes.push_back(neighbors[idx]);
                tree.edges.push_back(Edge{neighbors[idx],root});
            }
            std::cout << "motif created" << std::endl;
            printMotif(tree);
            // push subtree to queue if new submotifs can be generated
            if (tree.nodes.size() != MAX_NODES)
            {
                subtreeQueue.push(tree);  // push subtree to queue
            }

            // increment gdv[node][orbit] for every orbit found for each node in motif
            updateGDV(tree,orbitmatcher,gdv);  
            std::cout << "gdv updated" << std::endl;

            // generate next combination
            combo_gen.next();
        } 


        // main loop
        while (!subtreeQueue.empty())
        {   
            std::cout << subtreeQueue.size() << " trees in queue" << std::endl;

            // pop motif from queue
            Motif tree = subtreeQueue.front();
            std::cout << "Current tree: ";
            printMotif(tree);

            // get neighboring edges for each endnode
            neighboringEdges.clear();
            for (Node v: tree.endnodes)
            {
                for (Node w : g.neighbors(v))
                {
                    // we only want edges where w is higher than root and w isnt part of the motif already
                    if (w > root && std::find(tree.nodes.begin(),tree.nodes.end(),w) == tree.nodes.end())
                    {
                        neighboringEdges.push_back(Edge{w,v});
                    }
                }
            }
            std::cout << neighboringEdges.size() << " neighboring edges: [";
            for (auto edge: neighboringEdges)
                std::cout << "(" << edge.u << "," << edge.v << "), ";
            std::cout << "]" << std::endl;
            
            // skip this tree if it has no valid neighboring edges
            if (neighboringEdges.size() == 0)
            {
                subtreeQueue.pop();
                continue;
            }

            // setup combination generator for neighboring edges
            AllCombinationGenerator edge_combo_gen(neighboringEdges.size(), std::min((int) neighboringEdges.size(), (int) (MAX_NODES - tree.nodes.size())));

            // 1. create new subtree using tree and an edge combination
            // 2. update queue with new subtree
            // 3. update gdv using subtree
            // 4. update gdv using subtree's induced backedges
            while(!edge_combo_gen.done)
            {
                std::cout << "edge combo " << edge_combo_gen.combo_cnt << ": [";
                for (auto idx: edge_combo_gen.indices)
                    std::cout << "(" << neighboringEdges[idx].u << "," << neighboringEdges[idx].v << "), ";
                std::cout << "]" << std::endl;

                if (skip_combo(tree,edge_combo_gen.indices,neighboringEdges))
                {
                    std::cout << "skipping combo" << std::endl;
                    edge_combo_gen.next();
                    continue;
                }

                //create new subtree
                Motif newTree;
                newTree.root = tree.root;
                // copy edges from parent
                for (size_t i=0; i<tree.edges.size(); i++)
                {
                    newTree.edges.push_back(tree.edges[i]); 
                }
                // copy nodes from parent
                for (Node node: tree.nodes)
                {
                    newTree.nodes.push_back(node); 
                }

                // NOTE: should newTree include parent's endnodes ???
                // // copy endnodes from parent
                // for (Node node: tree.endnodes)
                // {
                //     newTree.endnodes.push_back(node); 
                // }

                // add new edges
                for (auto edge_idx : edge_combo_gen.indices)
                {
                    e = neighboringEdges[edge_idx];

                    newTree.edges.push_back(e);  // add new edge to motif
                    newTree.nodes.push_back(e.u);  // add new node to motif
                    newTree.endnodes.push_back(e.u); // mark new node as an endnode
                    
                } // end for adding new edges

                // add new tree to queue
                subtreeQueue.push(newTree);

                // update gdv for motifs generated with Tree and backedges
                update_gdv_backedges(newTree, g, gdv);
                
                // update the gdv using the new tree
                updateGDV(newTree,orbitmatcher,gdv);  

                std::cout << "gdv updated 2" << std::endl;

                // generate next combination
                edge_combo_gen.next();
            }  // end while new subtree
            subtreeQueue.pop();  // remove tree from queue and delete tree motif
        }  // end while mainloop
    } // end for
    return;
}

void printGDVTable(const std::vector<std::array<int,72>> gdv_table)
{
    std::cout << "final gdv" << std::endl;
    for (size_t i = 0; i < gdv_table.size(); i++)
    {
        std::cout << i << " : ";
        for (auto freq: gdv_table[i])
        {
            std::cout << freq << " ";
        }
        std::cout << std::endl;
    }
    return;
}

int main()
{
    
    //Edge edges[] = {{2,3}, {2,4}, {3,5}, {3,6}, {4,5}, {5,6}};
    Edge edges[] = {{0,1}, {0,2}, {1,3}, {1,4}, {2,3}, {3,4}};
    Edge petersen[] = {{0,1},{0,4},{0,5},{1,2},{1,6},{2,3},{2,7},{3,4},{3,8},{4,9},{5,7},{5,8},{6,8},{6,9},{7,9}};
    Edge karate[] = {
        {0,1},{0,2},{0,3},{0,4},{0,5},{0,6},{0,7},{0,8},{0,9},{0,10},{0,11},{0,12},{0,13},{0,14},{0,15},{0,16},{1,0},{1,14},{1,16},{1,18},{1,19},{2,0},{2,17},{2,18},{2,19},{2,20},{3,0},{3,16},{4,0},{4,16},{5,0},{5,16},{6,0},{6,17},{6,18},{7,0},{7,16},{8,0},{8,16},{9,0},{9,11},{9,13},{9,16},{9,32},{10,0},{10,13},{11,0},{11,9},{11,19},{11,33},{12,0},{12,15},{12,19},{13,0},{13,9},{13,10},{13,16},{14,0},{14,1},{14,16},{14,17},{15,0},{15,12},{15,16},{15,18},{15,32},{15,33},{16,0},{16,1},{16,3},{16,4},{16,5},{16,7},{16,8},{16,9},{16,13},{16,14},{16,15},{16,19},{17,2},{17,6},{17,14},{17,18},{17,19},{17,20},{17,24},{17,30},{17,31},{18,1},{18,2},{18,6},{18,15},{18,17},{18,19},{18,20},{18,21},{18,22},{18,23},{18,24},{18,26},{18,27},{18,28},{18,30},{18,31},{19,1},{19,2},{19,11},{19,12},{19,16},{19,17},{19,18},{19,20},{19,24},{19,25},{20,2},{20,17},{20,18},{20,19},{20,24},{20,28},{21,18},{21,23},{21,26},{22,18},{22,23},{22,26},{22,29},{23,18},{23,21},{23,22},{23,29},{24,17},{24,18},{24,19},{24,20},{25,19},{26,18},{26,21},{26,22},{27,18},{28,18},{28,20},{29,22},{29,23},{30,17},{30,18},{31,17},{31,18},{32,9},{32,15},{32,33},{33,11},{33,15},{33,32}
    };
    // Graph graph(edges,5,6);
    Graph graph(petersen,10,15);
    //Graph graph(karate,34,154);
    std::cout << "Graph: ";
    graph.printGraph();
    std::cout << "creating orbitmatcher" << std::endl;
    OrbitMatcher orbitmatcher;
    std::cout << "finished creating orbitmatcher" << std::endl;
    // read_orbit_table(orbitmatcher);
    std::vector<std::array<int,72>> gdvs;
    for (int i = 0; i < 5; i++)
    {
        gdvs.push_back(std::array<int,72>());
        gdvs[i].fill(0);
    }
    std::cout << "starting computation" << std::endl;
    computeGDV(graph,orbitmatcher,gdvs);
    std::cout << "ending computation" << std::endl;
    printGDVTable(gdvs);
}