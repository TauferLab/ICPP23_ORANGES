#ifndef GPUBFS_HPP
#define GPUBFS_HPP

#include "raw_vecs.hpp"
#include "structure_defs.hpp"
#include <iostream>

using namespace std;

// DFS algorithm for A_Network that returns list of visited nodes.
FIDO_HOST_DEVICE void a_network_dir_dfs_raw(A_Network_raw& network, intvec& visited)
{
    // Set up BFS
    // visited.clear();
    clear_intvec(visited);
    // Max size 6
    // vector<int> queue;
    intvec queue = new_intvec(6);
    // queue.push_back(0);
    pushback_intvec(queue, 0);
    int head;

    // Make variables for containing checks
    bool bfs_visited;
    int visited_counter;
    bool is_neighbor;
    int neighbor_counter;

    // Loop through rest of network
    // while (!queue.empty()) {
    while (queue.veclen > 0)
    {
        // head = queue.front();
        // queue.erase(queue.begin());
        // head = queue.back();
        head = queue.vec[queue.veclen - 1];
        // queue.pop_back();
        popback_intvec(queue);
        int found_neighbor = -1;

        // Search through possible neighbors of the head
        // for (int i = 0; i < network.size(); i++) {
        for (int i = 0; i < network.nodes_len; i++)
        {
            // Check if network[i] corresponds to an in or out neighbor of network[head]
            is_neighbor      = false;
            neighbor_counter = 0;
            // while (!is_neighbor && (neighbor_counter < network[i].ListW.size())) {
            while (!is_neighbor && (neighbor_counter < network.vec[i].ListW.veclen))
            {
                // if (network[i].ListW[neighbor_counter].first == network[head].Row) {
                if (network.vec[i].ListW.vec[neighbor_counter].first == network.vec[head].Row)
                {
                    is_neighbor = true;
                    // found_neighbor = i;
                }
                neighbor_counter += 1;
            }
            neighbor_counter = 0;
            // while (!is_neighbor && (neighbor_counter < network[head].ListW.size())) {
            while (!is_neighbor && (neighbor_counter < network.vec[head].ListW.veclen))
            {
                if (network[head].ListW[neighbor_counter].first == network[i].Row)
                {
                    is_neighbor = true;
                    // found_neighbor = i;
                }
                neighbor_counter += 1;
            }


            // if (is_neighbor) {
            //	cout << network[i].Row << " is a neighbor of " << network[head].Row << endl;
            //} else {
            //	cout << network[i].Row << " is not a neighbor of " << network[head].Row << endl;
            //}

            if (is_neighbor)
            {
                // Determine if neighbor is already visited
                bfs_visited     = false;
                visited_counter = 0;
                // while (!bfs_visited && (visited_counter < visited.size())) {
                while (!bfs_visited && (visited_counter < visited.veclen))
                {
                    if (network.vec[i].Row == visited.vec[visited_counter])
                    {
                        bfs_visited = true;
                    }
                    visited_counter += 1;
                }

                // If neighbor is not already visited, add neighbor to queue
                if (!bfs_visited)
                {
                    // queue.push_back(i);
                    pushback_intvec(queue, i);
                }
            }
        }
        //  }
        // visited.push_back(network[head].Row);
        pushback_intvec(visited, network.vec[head].Row);
    }
    delete_intvec(queue);
    return;
}


// Returns the shortest path length to node from every node in network. Path lengths are stored in
// distance and indexed by the index of network.
FIDO_HOST_DEVICE void dir_dfs_shortest_paths_raw(int node, A_Network_raw& network, intvec& distance)
{
    // Set up distance vector
    // distance.clear();
    int node_index;
    // distance.resize(network.size());
    distance.veclen = network.nodes_len;
    for (int i = 0; i < distance.veclen; i++)
    {
        distance.vec[i] = distance.veclen + 2;
        if (network.vec[i].Row == node)
        {
            node_index = i;
        }
    }
    distance.vec[node_index] = 0;

    // Declare needed variables for algorithm
    // vector<int> queue;
    intvec queue = new_intvec(6);
    // queue.push_back(node_index);
    pushback_intvec(queue, node_index);
    int current_node;
    int neighbor_index;

    bool is_neighbor;
    int neighbor_counter;

    // Loop through rest of network
    // while (!queue.empty()) {
    while (queue.veclen > 0)
    {
        // current_node = queue.front();
        // queue.erase(queue.begin());
        // current_node = queue.back();
        current_node = queue.vec[queue.veclen - 1];
        // queue.pop_back();
        popback_intvec(queue);

        // Update system if path to neighbor is longer than the distance to the current node + 1
        // for (int i = 0; i < network.size(); i++) {
        for (int i = 0; i < network.nodes_len; i++)
        {
            // Check if network[i] corresponds to an in or out neighbor of network[head]
            is_neighbor      = false;
            neighbor_counter = 0;
            // while (!is_neighbor && (neighbor_counter < network[i].ListW.size())) {
            while (!is_neighbor && (neighbor_counter < network.vec[i].ListW.veclen))
            {
                // if (network[i].ListW[neighbor_counter].first == network[current_node].Row) {
                if (network.vec[i].ListW.vec[neighbor_counter].first ==
                    network.vec[current_node].Row)
                {
                    is_neighbor = true;
                }
                neighbor_counter += 1;
            }
            neighbor_counter = 0;
            // while (!is_neighbor && (neighbor_counter < network[current_node].ListW.size())) {
            while (!is_neighbor && (neighbor_counter < network.vec[current_node].ListW.veclen))
            {
                // if (network[current_node].ListW[neighbor_counter].first == network[i].Row) {
                if (network.vec[current_node].ListW.vec[neighbor_counter].first ==
                    network.vec[i].Row)
                {
                    is_neighbor = true;
                }
                neighbor_counter += 1;
            }

            // if ((distance[i] > distance[current_node] + 1) && is_neighbor) {
            if ((distance.vec[i] > distance.vec[current_node] + 1) && is_neighbor)
            {
                distance.vec[i] = distance.vec[current_node] + 1;
                // queue.push_back(i);
                pushback_intvec(queue, i);
            }
        }
    }
    delete_intvec(queue);
    return;
}

#endif
