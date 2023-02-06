
#ifndef ESSENS_GRAPH_
#define ESSENS_GRAPH_

#include<vector>
#include<map>
#include<set>
#include<unordered_set>
#include<queue>
#include<iostream>
#include<memory>

// TODO: use smart pointers

namespace ESSENS
{

typedef unsigned int Node;

struct Edge
{
    Node u;
    Node v;
};

bool operator== (Edge first, Edge second)
{
    return std::make_pair(first.u,first.v) < std::make_pair(second.u,second.v);
}

bool operator< (Edge first, Edge second) 
{ 
    return std::make_pair(first.u,first.v) < std::make_pair(second.u,second.v);
}

class Motif
{
    public:
    Node root; // Node is basically int. Edge is basically pair<Node,Node>
    //std::vector<std::pair<Node,bool>> nodes; // pair with first int the node id and the second whether the node is a leaf
    //std::map<Node,bool> nodes;  // nodes and whether the node is a leaf
    std::vector<Node> nodes;
    std::vector<Node> endnodes; 
    //std::array<Node,5> Nodes;
    //std::array<bool, 5> endNodes;
    //std::vector<std::shared_ptr<Edge>> edges;
    std::vector<Edge> edges;
    //std::vector<Edge*> edges;
    std::vector<Node> get_neighbors(Node node) const;
    std::map<Node,int> get_degrees() const;
    Motif() {}
    Motif(  int root_, 
            std::vector<Node> nodes_, 
            std::vector<Node> endnodes_, 
            std::vector<Edge> edges_)
            :root(root_),nodes(nodes_),endnodes(endnodes_),edges(edges_){};
};


// Motif::Motif(Graph g, Motif parent, std::vector<std::shared_ptr<Edge>> new_edges)
// {
    
// }

// get neighbors of the node in the motif
std::vector<Node> Motif::get_neighbors(Node node) const
{
    // iterate through edges in motif
    std::vector<Node> neighbors;
    for (Edge edge: edges)
    {
        // if node is part of an edge, we add the other node from the edge to the neighbors vector
        if (edge.u == node)
        {
            neighbors.push_back(edge.v);
        } else if (edge.v == node)
        {
            neighbors.push_back(edge.u);
        }
    } // end for
    // std::vector<std::shared_ptr<Edge>>::const_iterator it;
    // for (it = edges.begin(); it != edges.end(); it++)
    // {
    //     // if node is part of an edge, we add the other node from the edge to the neighbors vector
    //     if ((*it)->u == node)
    //     {
    //         neighbors.push_back((*it)->v);
    //     } else if ((*it)->v == node)
    //     {
    //         neighbors.push_back((*it)->u);
    //     }
    // } // end for
    return neighbors;
}

// get degree of 
std::map<Node,int> Motif::get_degrees() const
{
    std::map<Node,int> degrees;
    for (Node node: nodes)
    {
        degrees.insert({node,0});
    }
    // iterate through edges in motif
    for (Edge edge: edges)
    {
        //degrees.at(edge.u)++;
        degrees[edge.u]++;
        degrees[edge.v]++;
    }
    // std::vector<std::shared_ptr<Edge>>::const_iterator it;
    // for (it = edges.begin(); it != edges.end(); it++)
    // {
    //     // increment 
    //     degrees.at((*it)->u)++;
    //     degrees.at((*it)->v)++;
    // } // end for
    return degrees;
}
// class Edge
// {
//     public:
//         Node u;
//         Node v;

//         Edge(Node a, Node b)
//         {
//             this->u = a;
//             this->v = b;
//         }
// };

class Graph
{
    public:
    int N;
    int M;
    std::set<Node> Nodes;
    std::map<Node,std::vector<Node>> adjList;
    //std::vector<Edge> edgelist;

    //Graph();
    Graph(const Edge edges[], int n, int m);
    void addEdge(const Edge &e);
    inline std::vector<Node> neighbors(const Node u) const;
    void edgelist(std::vector<Edge> &edges) const;
    // void inducedSubgraph(const std::vector<Node> &nodes, Graph &subgraph);
    // void traverse(Graph g, Node root, std::vector<Node> &traversal, std::map<Node,Node> &parents);
    void printGraph() const;
    // void all_pair_shortest_path(std::vector<std::vector<std::vector<Node>>> paths);
    // void single_source_shortest_path(Node node, std::vector<std::vector<Node>> paths);
};


Graph::Graph(const Edge edges[], int n, int m):N(n),M(m)
{
    int i;
    for (i = 0; i < M; i++)
    {
        addEdge(edges[i]);
        // std::cout << "added edge (" << edges[i].u << ", " << edges[i].v << ")" << std::endl;
    }
}

void Graph::addEdge(const Edge &edge)
{
    // insert empty vector at u in map if it doesnt already exist
    std::pair<std::map<Node, std::vector<Node>>::iterator,bool> ret;
    ret = adjList.insert(std::map<Node,std::vector<Node>>::value_type(edge.u, std::vector<Node>()));
    ret.first->second.push_back(edge.v);  // add v to u's adj list
    Nodes.insert(edge.u);  // add u to the set of nodes

    // insert empty vector at v in map if it doesnt already exist
    ret = adjList.insert(std::map<Node,std::vector<Node>>::value_type(edge.v, std::vector<Node>()));
    ret.first->second.push_back(edge.u);  // add u to v's adj list
    Nodes.insert(edge.v);  // add v to the set of nodes
}

inline std::vector<Node> Graph::neighbors(Node node) const
{
    return adjList.find(node)->second;
}

void Graph::edgelist(std::vector<Edge> &edges) const
{
    Node u;
    Edge edge;
    std::set<Edge> seen;
    for (auto const &node_neighbors: adjList)
    {
        u = node_neighbors.first;
        for (Node v: node_neighbors.second)
        {
            edge = Edge{u,v};
            if (seen.find(edge) == seen.end())
            {
                edges.push_back(edge);
                seen.insert(edge);
                seen.insert(Edge{v,u});
            }
        }
    }
}

// void Graph::traverse(Graph g, Node root, std::vector<Node> &traversal, std::map<Node,Node> &parents)
// {
//     //bfs
//     std::unordered_set<Node> seen;
//     std::queue<Node> q;
//     Node curr, neighbor;

//     q.push(root);
//     seen.insert(root);

//     while(!q.empty())
//     {
//         curr = q.front();
//         traversal.push_back(curr);
//         q.pop();
//         for(Node neighbor: g.neighbors(curr))
//         {
//             if (seen.find(neighbor) == seen.end())
//             {
//                 parents.insert({curr, neighbor});
//                 q.push(neighbor);
//             }
//         }
//     }
// }

/*
void Graph::inducedSubgraph(Graph g, const std::vector<Node> &nodes, Graph &subgraph)
{
    std::unordered_set<Node> seen;
    int i,j;
    Node u;

    for (i = 0; i < nodes.size(); i++)
    {
        u = nodes[i];
        for (Node v: g.neighbors(u))
        {
            if (seen.find(v) != seen.end())
            {
                // prevent adding (v,u) if we already have added (u,v)
                continue;
            }
            for (j = i; j < nodes.size(); j++)
            {
                if (nodes[j] == v) 
                {
                    subgraph.addEdge(Edge{u,v});
                }
            }
        }
    }
}
*/

//void all_source_shortest_lengths(Motif motif, std::vector<std::vector<std::vector<Node>>> paths)
void all_pair_shortest_lengths(const Motif &motif, std::map<std::pair<Node,Node>,int> &distances)
{
    // bfs
    std::unordered_set<Node> seen;  // visited nodes
    std::queue<std::pair<Node,int>> q;  // queue of (node,distance)
    Node curr;
    std::pair<Node, int> node_dist;  // pair of node and its distance to root
    // std::vector<Node> neighbors;
    int distance = 0;

    // iterate through nodes and use them as the source of a bfs
    for (Node root: motif.nodes)
    {
        // std::cout << std::endl << "root=" << root << std::endl;
        // bfs
        // initialize set and queue
        seen.clear();
        q.push({root,0});  
        seen.insert(root); 

        while(!q.empty())
        {
            // pop (node,distance) from queue
            node_dist = q.front();
            q.pop();

            curr = node_dist.first;  // current node
            distance = node_dist.second + 1;  // new distance
            // std::cout << "(" << curr << ", d=" << node_dist.second << ")" << std::endl;
            // // get neighbors
            // neighbors.clear();
            // motif.get_neighbors(curr, neighbors);

            // iterate over neighbors and update queue accordingly
            for (Node neighbor: motif.get_neighbors(curr))
            {
                // std::cout << "curr=" << curr << " neigbhor=" << neighbor << std::endl;
                // add to queue if neighbor not yet seen
                if (seen.find(neighbor) == seen.end())
                {
                    q.push({neighbor,distance});  // add (neighbor,distance) to queue
                    distances.at({root,neighbor}) = distance; // update distances[(curr,neighbor)]
                    seen.insert(neighbor);  // add neighbor to the set of visited nodes
                }
            }  // end for neighbors 
        }  // end while
    } //end for
    return;
}

void kcore_counts(const std::map<Node,int> &count)
{
    // iterative dfs to get kcores
}

void Graph::printGraph() const
{
    std::vector<Edge> edges;
    edgelist(edges);
    for(Edge e: edges)
    {
        std::cout << "(" << e.u << "," << e.v << ") ";
    }
    std::cout << std::endl;
}

void printMotif(const Motif &motif)
{
    std::cout << "Root=" << motif.root << " Nodes= [";
    for (Node node: motif.nodes)
    {
        std::cout << node << ",";
    }
    std::cout << "] Endnodes= [";
    for (Node endnode: motif.endnodes)
    {
        std::cout << endnode << ",";
    }
    std::cout << "] Edges= [";
    for (Edge edge: motif.edges)
    {
        std::cout << "(" << edge.u << "," << edge.v << "),";
    }
    std::cout << "]" << std::endl;
}

} // essens namespace
#endif