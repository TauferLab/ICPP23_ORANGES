#ifndef ORBIT_MATCH_
#define ORBIT_MATCH_

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <array>
#include <stdexcept>
#include <algorithm>

#include "graph.hpp"

using namespace ESSENS;

class OrbitMatcher
{
    public:
        struct OrbitInfo
        {
            int orbit;
            int graphlet;
            std::array<int, 5> degree_vector;
            std::array<int, 5> sp_vector;
        };

        std::array<OrbitInfo, 73> OrbitTable;  // table of orbits and their features

        std::vector<int> get_orbits(const Motif &motif) const;

        OrbitMatcher()
        {
            for (auto it = OrbitTable.begin(); it != OrbitTable.end(); it++)
            {
                *it = {0,0,std::array<int,5>(),std::array<int,5>()};
                (*it).degree_vector.fill(0);
                (*it).sp_vector.fill(0);
            }
            read_orbit_info("orbitinfo.tsv");
        }
    
    private:
        void read_orbit_info(const std::string fname);
        int match_orbit(std::array<int,5> deg_sig, std::array<int,5> sp_sig) const;
};

// read tsv file into OrbitTable
void OrbitMatcher::read_orbit_info(const std::string fname)
{
    int i = 0;
    int j = 0;
    std::ifstream fin(fname);
    std::string row, col, tmp;

    std::getline(fin, row); // skip the 1st line
    // iterate over the lines in the file
    while(getline(fin, row))
    {
        std::istringstream ss(row);

        // Graphlet
        std::getline(ss, col, '\t');
        OrbitTable[i].graphlet = stoi(col);

        // Orbit
        std::getline(ss, col, '\t');
        OrbitTable[i].orbit = stoi(col);

        // Degree Vector
        std::getline(ss, col, '\t');
        col = col.substr(1,col.size()-2);  // remove '[' and ']' from string
        std::istringstream ss2(col);
        j = 0;
        // Iterate through comma seperated integers and fill array
        while(getline(ss2, tmp, ','))
        {
            OrbitTable[i].degree_vector[j] = stoi(tmp);
            ++j;
        }
        // sort the deg vector so the zeros go in the beginning if len(deg_vec) < 5
        std::sort(OrbitTable[i].degree_vector.begin(), OrbitTable[i].degree_vector.end());
        
        std::getline(ss, col, '\t');
        col = col.substr(1,col.size()-2);  // remove '[' and ']' from string
        std::istringstream ss3(col);
        j = 0;
        // Iterate through comma seperated integers and fill array
        while(getline(ss3, tmp, ','))
        {
            OrbitTable[i].sp_vector[j] = stoi(tmp);
            ++j;
        }

        ++i;
    }  //endwhile
    return;
}  //end function read_orbit_info()

// given degree and shortest path signatures, return the corresponding orbit from the OrbitTable
int OrbitMatcher::match_orbit(std::array<int,5> deg_sig, std::array<int,5> sp_sig) const
{
    // iterate through the orbit table
    for (OrbitInfo orbit_info: OrbitTable)
    {
        // return when deg_sig and sp_sig match to an orbit's degree and shortest path vectors
        if (orbit_info.degree_vector == deg_sig && orbit_info.sp_vector == sp_sig)
        {
            return orbit_info.orbit;
        }
    }
    return -1;  // -1 if no orbit found
}

// given motif, generate the shortest path signatures for each node in motif as root
void get_sp_sigs(const Motif &motif, std::map<Node,std::array<int,5>> &sp_sigs)
{
    // get total number of shortest paths from root of distance x
    // performs bfs on every node as source
    std::map<Node,int> seen;  // visited nodes and their level
    std::queue<std::pair<Node,int>> q;  // queue of (node,distance)
    Node curr;
    std::pair<Node, int> node_dist;  // pair of node and its distance to root
    int distance = 0;

    // iterate through nodes and use them as the source of a bfs
    for (Node root: motif.nodes)
    {
        // bfs
        // initialize set and queue
        seen.clear();
        q.push({root,0});  // root with distance 0
        seen.insert({root,0});

        while(!q.empty())
        {
            // get front of queue 
            node_dist = q.front();  // (node,distance)
            curr = node_dist.first;  // current node
            distance = node_dist.second + 1;  // new distance

            // iterate over neighbors and update queue accordingly
            for (Node neighbor: motif.get_neighbors(curr))
            {
                // add to queue if neighbor not yet seen
                if (seen.find(neighbor) == seen.end())
                {
                    q.push({neighbor,distance});  // add (neighbor,distance) to queue
                    sp_sigs[root][distance]++; // increment the shortest path vector
                    seen.insert({neighbor,distance});  // add (neighbor,distance) to the visited set
                } else if (distance == seen.find(neighbor)->second)  // if neighbor has been seen and its distance is the same as current
                {
                    // this is another shortest path, so we increment the shortest path counter
                    sp_sigs[root][distance]++;
                }
            }  // end for neighbors 
            q.pop();  // remove (node,distance) from the queue
        }  // end while
    } //end for

    return;
}

// finds the correspoding orbit for each node in the motif and increments gdv accordingly
//void OrbitMatcher::get_orbits(const Motif &motif, std::map<Node,std::array<int,5>> &gdv) const
//void OrbitMatcher::update_gdv(const Motif &motif, std::vector<std::array<int,72>> &gdv) const
std::vector<int> OrbitMatcher::get_orbits(const Motif &motif) const
{
    // initialize variables
    std::array<int, 5> deg_sig;  // degree signature of motif
    std::map<Node, std::array<int,5>> sp_sigs;  // map of {node: shortest_path_vector}
    int orbit;
    int N = motif.nodes.size();
    std::vector<int> orbits;
    
    // std::cout << "Motif: ";
    // printMotif(motif);

    if(N > 5)
    {
        throw std::invalid_argument("Motif has more than 5 nodes");
    }
    
    deg_sig.fill(0);  // initialize deg_sig
    // initialize sp_sigs
    for (Node node: motif.nodes)
    {
        sp_sigs.insert({node, std::array<int,5>()});
        sp_sigs[node].fill(0);
        sp_sigs[node][0] = 1;
    }

    // std::cout << "spsig " << std::endl;;
    // for (auto node_sp_sig: sp_sigs)
    // {
    //     std::cout << "(" << node_sp_sig.first << "): ";
    //     for (auto freq: node_sp_sig.second)
    //     {
    //         std::cout << freq << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // calculate degree signature
    size_t i = 0;
    for (auto const& node_degree: motif.get_degrees())
    {
        deg_sig[i] = node_degree.second;  // store the degree of the node at index i
        i++;
    }
    std::sort(deg_sig.begin(), deg_sig.end());

    // calculate sp signature. sp[x] = number of shortest paths of length x
    get_sp_sigs(motif, sp_sigs);

    // match deg_sig and sp_sig to orbit
    for (Node node: motif.nodes)
    {
        orbit = match_orbit(deg_sig,sp_sigs[node]);
        // std::cout << "Node=" << node << " Orbit=" << orbit;
        // std::cout << " degsig=[";
        // for (auto deg: deg_sig)
        // {
        //     std::cout << deg << ",";
        // }
        // std::cout << "] spsig= [";
        // for (auto freq: sp_sigs[node])
        // {
        //     std::cout << freq << ",";
        // }
        // std::cout << "]" << std::endl;
        orbits.push_back(orbit);
    }
    return orbits;
}


void read_orbit_table(OrbitMatcher orbit_matcher)
{
    auto orbit_table = orbit_matcher.OrbitTable;
    for (auto orbit_info: orbit_table)
    {
        std::cout << orbit_info.orbit << " " << orbit_info.graphlet << " [";
        for (auto deg: orbit_info.degree_vector)
        {
            std::cout << deg << ",";
        }
        std::cout << "] [";
        for (auto sp: orbit_info.sp_vector)
        {
            std::cout << sp << ",";
        }
        std::cout << "]" << std::endl;
    }
}

#endif