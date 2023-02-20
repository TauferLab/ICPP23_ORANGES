
#ifndef ESSENS_GRAPH_
#define ESSENS_GRAPH_

#include<vector>
#include<map>
#include<set>
#include<unordered_set>
#include<unordered_map>
#include<queue>
#include<iostream>
#include<memory>
#include <string>
#include <fstream>
#include <sstream>

// TODO: use smart pointers

namespace ESSENS
{

typedef unsigned int Node;

struct Edge
{
    Node u;
    Node v;
    ~Edge() {}
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
    std::vector<Node> nodes;
    std::vector<Node> endnodes; 
    // std::array<Node,5> nodes;
    // std::array<Node,5> endnodes;
    // size_t num_nodes;
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
            :root(root_),nodes(nodes_),endnodes(endnodes_),edges(edges_)
    {}
};

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
    return degrees;
}

class Graph
{
    public:
    int N;
    int M;
    std::set<Node> Nodes;
    std::map<Node,std::vector<Node>> adjList;

    Graph(const std::vector<Edge> &edges);
    Graph(const Edge edges[], int n, int m);
    void addEdge(const Edge &e);
    inline std::vector<Node> neighbors(const Node u) const;
    void edgelist(std::vector<Edge> &edges) const;
    // void inducedSubgraph(const std::vector<Node> &nodes, Graph &subgraph);
    // void traverse(Graph g, Node root, std::vector<Node> &traversal, std::map<Node,Node> &parents);
    void printGraph() const;
    // void all_pair_shortest_path(std::vector<std::vector<std::vector<Node>>> paths);
    // void single_source_shortest_path(Node node, std::vector<std::vector<Node>> paths);
    void relabel_graph(std::vector<Node> const &ordering);
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

Graph::Graph(const std::vector<Edge> &edges)
{
    // assumes undirected and unique
    for (Edge edge: edges)
    {
        addEdge(edge);
    }
    N = Nodes.size();
    M = edges.size();
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

void kcore_ordering(Graph graph, std::vector<Node> &ordered)
{
    // given graph, orders the nodes into a vector based on their kcore number
    // aka degeneracy ordering
    // Matula & Beck (1983)
    std::unordered_map<Node, int> degrees;  // the number of neighbors of v that are not already in output
    std::vector<std::vector<Node>> buckets;  // degree-indexed bucket queue
    int max_deg = -1;  // maximum degree of the graph
    for (auto const &node_pair: graph.adjList)
    {
        int deg = node_pair.second.size();
        degrees[node_pair.first] = deg;  // initialize to degree of the node
        if (max_deg < deg)
        {
            max_deg = deg;
        }
    }

    for (int i = 0; i <= max_deg; i++)
    {
        buckets.emplace_back();
    }
    // initialize buckets for which deg(v)=i
    for (auto const &node_pair: graph.adjList)
    {
        int deg = node_pair.second.size();
        buckets[deg].push_back(node_pair.first);
    }
    
    size_t k = 0;

    // iterate N times
    for (int x = 0; x < graph.N; x++)
    {
        // find nonempty degrees[i]
        size_t i;
        for (i = 0; i < buckets.size() && buckets[i].size() <= 0; i++);

        k = std::max(k,i);
        // pop node from bucket[i] and add it to ordered
        Node u = buckets[i].back();
        ordered.push_back(u);
        buckets[i].pop_back();

        // update neighbors
        for (Node neighbor: graph.adjList[u])
        {
            // check if neighbor isnt already in ordered
            if ( std::find(ordered.begin(),ordered.end(),neighbor) == ordered.end())
            {
                int deg = degrees[neighbor];  // current degree of neighbor
                degrees[neighbor]--;  // decrement degree of neighbor
                // move neighbor from current bucket to new bucket
                buckets[deg-1].push_back(neighbor);  // add neighbor to new bucket
                auto it = std::find(buckets[deg].begin(),buckets[deg].end(),neighbor);
                buckets[deg].erase(it);  // remove neighbor from its current bucket
            }  // end if
        }  // end for updating neighbors
    } // end for
}  // end function core numbers


void Graph::relabel_graph(std::vector<Node> const &ordering)
{
    // relabels the nodes in the graph where ordering[i] is the old and i is the new label

    // initialize map for faster access
    std::unordered_map<Node,Node> node_ordering; // {old label : new label}
    for (std::size_t i = 0; i != ordering.size(); i++)
    {
        node_ordering[ordering[i]] = (Node) i;
    }

    std::map<Node,std::vector<Node>> new_adjList;  // new adjacency list
    for (auto node_pair: node_ordering)
    {
        Node u = node_pair.second;
        Node old_u = node_pair.first;
        std::vector<Node> old_neighbors = adjList[old_u];
        std::vector<Node> new_neighbors;
        for (Node old_v: old_neighbors)
        {
            Node new_v = ordering.at(old_v);
            new_neighbors.push_back(new_v);
        }
        new_adjList[u] = new_neighbors;
    }
    adjList = new_adjList;
}  // end function relabel_graph


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

std::vector<Edge> read_graph(std::string fname)
{
    // assumes first line has number of nodes, and number of edges
    // every subsequent line is an edge
    // note that first node should be labelled 0
    std::vector<Edge> edgelist;
    std::ifstream ifs(fname);
    if (!ifs.is_open())
        std::cout << fname << " is not open!" << std::endl;

    std::string line, word1, word2;

    getline(ifs, line); // skip first line

    while(getline(ifs, line))
    {
        std::stringstream str(line);
        std::getline(str,word1,' ');
        std::getline(str,word2,' ');
        edgelist.push_back({(Node)stoi(word1),(Node)stoi(word2)});
    }
    return edgelist;
}


} // essens namespace
#endif