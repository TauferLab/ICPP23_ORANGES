import networkx as nx
import time
import concurrent.futures
import os
import numpy as np
from xml.dom import minidom

start = time.perf_counter()

class OrbitMetric:
    def __init__(self,orbit_number, orbit_degree , orbit_distance):
        self.orbit_degree = orbit_degree
        self.orbit_distance = orbit_distance
        self.orbit_number = orbit_number

class GDVMetric:
    def __init__(self,GDV,node):
        self.GDV = GDV
        self.node = node

# This function helps to create a graph from text file with edge list with nodes seperated by space
def graph_creation_from_file_txt(path):
    G = nx.Graph()
    file1 = open(path, 'r') 
    Lines = file1.readlines() 
    for line in Lines:
        edge = line.strip().split(' ')
        G.add_edge(str(edge[0]),str(edge[1]))
    return G

def orbit_list_creation(n_count = None):
    # Total 72 orbits with 0 to 71 indices each with 5 distance nodes (0 to 5 distance)
    # Currently only 21 orbits are considered.
    orbits = []
    # If node count is not specified, then return all orbits
    if n_count is None:
        # orbits with 2 nodes
        orbits.append(OrbitMetric(0,[1,1],[1,1,0,0,0,0]))
        # orbits with 3 nodes
        orbits.append(OrbitMetric(1,[1,2,1],[1,1,1,0,0,0]))
        orbits.append(OrbitMetric(2,[1,2,1],[1,2,0,0,0,0]))
        orbits.append(OrbitMetric(3,[2,2,2],[1,2,0,0,0,0]))
        # orbits with 4 nodes
        orbits.append(OrbitMetric(4,[1,2,2,1],[1,1,1,1,0,0]))
        orbits.append(OrbitMetric(5,[1,2,2,1],[1,2,1,0,0,0]))
        orbits.append(OrbitMetric(6,[1,1,1,3],[1,1,2,0,0,0]))
        orbits.append(OrbitMetric(7,[1,1,1,3],[1,3,0,0,0,0]))
        orbits.append(OrbitMetric(8,[2,2,2,2],[1,2,2,2,0,0]))
        orbits.append(OrbitMetric(9,[1,3,2,2],[1,1,2,2,0,0]))
        orbits.append(OrbitMetric(10,[1,3,2,2],[1,2,3,1,0,0]))
        orbits.append(OrbitMetric(11,[1,3,2,2],[1,3,2,0,0,0]))
        orbits.append(OrbitMetric(12,[2,2,3,3],[1,2,3,1,0,0]))
        orbits.append(OrbitMetric(13,[2,2,3,3],[1,3,1,2,0,0]))
        orbits.append(OrbitMetric(14,[3,3,3,3],[1,3,3,2,0,0]))
        # orbits with 5 nodes
        orbits.append(OrbitMetric(15,[1,2,2,2,1],[1,1,1,1,1,0]))
        orbits.append(OrbitMetric(16,[1,2,2,2,1],[1,2,1,1,0,0]))
        orbits.append(OrbitMetric(17,[1,2,2,2,1],[1,2,2,0,0,0]))
        orbits.append(OrbitMetric(18,[1,2,3,1,1],[1,1,1,2,0,0]))
        orbits.append(OrbitMetric(19,[1,2,3,1,1],[1,1,2,1,0,0]))
        orbits.append(OrbitMetric(20,[1,2,3,1,1],[1,2,2,0,0,0]))
        orbits.append(OrbitMetric(21,[1,2,3,1,1],[1,3,1,0,0,0]))
    elif n_count ==2:
        orbits.append(OrbitMetric(0,[1,1],[1,1,0,0,0,0]))
    elif n_count == 3:
        orbits.append(OrbitMetric(1,[1,2,1],[1,1,1,0,0,0]))
        orbits.append(OrbitMetric(2,[1,2,1],[1,2,0,0,0,0]))
        orbits.append(OrbitMetric(3,[2,2,2],[1,2,0,0,0,0]))
    elif n_count == 4:
        orbits.append(OrbitMetric(4,[1,2,2,1],[1,1,1,1,0,0]))
        orbits.append(OrbitMetric(5,[1,2,2,1],[1,2,1,0,0,0]))
        orbits.append(OrbitMetric(6,[1,1,1,3],[1,1,2,0,0,0]))
        orbits.append(OrbitMetric(7,[1,1,1,3],[1,3,0,0,0,0]))
        orbits.append(OrbitMetric(8,[2,2,2,2],[1,2,2,2,0,0]))
        orbits.append(OrbitMetric(9,[1,3,2,2],[1,1,2,2,0,0]))
        orbits.append(OrbitMetric(10,[1,3,2,2],[1,2,3,1,0,0]))
        orbits.append(OrbitMetric(11,[1,3,2,2],[1,3,2,0,0,0]))
        orbits.append(OrbitMetric(12,[2,2,3,3],[1,2,3,1,0,0]))
        orbits.append(OrbitMetric(13,[2,2,3,3],[1,3,1,2,0,0]))
        orbits.append(OrbitMetric(14,[3,3,3,3],[1,3,3,2,0,0]))
    elif n_count == 5:
        orbits.append(OrbitMetric(15,[1,2,2,2,1],[1,1,1,1,1,0]))
        orbits.append(OrbitMetric(16,[1,2,2,2,1],[1,2,1,1,0,0]))
        orbits.append(OrbitMetric(17,[1,2,2,2,1],[1,2,2,0,0,0]))
        orbits.append(OrbitMetric(18,[1,2,3,1,1],[1,1,1,2,0,0]))
        orbits.append(OrbitMetric(19,[1,2,3,1,1],[1,1,2,1,0,0]))
        orbits.append(OrbitMetric(20,[1,2,3,1,1],[1,2,2,0,0,0]))
        orbits.append(OrbitMetric(21,[1,2,3,1,1],[1,3,1,0,0,0]))

    return orbits

# Global variable created as orbit list is called most of the time
orbit_list = orbit_list_creation()


def GDV_calculation(vertex,Graph):
    GDV= [0]*len(orbit_list)
    samp = nx.single_source_shortest_path_length(Graph, vertex, cutoff=4)
    neighbour_list = list(samp.keys())
    neighbour_list.remove(vertex)
    from itertools import combinations 
    for i in range(1,5):
        comb = combinations(neighbour_list, i) 
        list_combinations = list(comb)
        for a in list_combinations:
            tuple_list = list(a)
            tuple_list.append(vertex)
            sub_graph = Graph.subgraph(tuple_list)
            if(nx.is_connected(sub_graph)):            
                new_graph = nx.Graph()  
                new_graph.add_edges_from(sub_graph.edges())
                new_graph_spl = nx.single_source_shortest_path_length(new_graph, vertex, cutoff=4)
                sub_graph_degree_list = new_graph.degree(new_graph.nodes())
                sub_graph_degree_list = list(sub_graph_degree_list)
                # Using sum() + values() 
                # Selective key values in dictionary 
                distance_list = [0]*6
                degree_list = [0]*len(new_graph.nodes())
                for g in range(0,len(degree_list)):
                    degree_list[g]= sub_graph_degree_list[g][1]
                for d in range(0,6):
                    res = sum(x == d for x in new_graph_spl.values())
                    distance_list[d]=res
                orbit_list_filter = orbit_list_creation(i+1)
                for n in range(0,len(orbit_list_filter)):
                    if(orbit_list_filter[n].orbit_distance==distance_list and sorted(list(orbit_list_filter[n].orbit_degree)) == sorted(degree_list) ):
                        GDV[orbit_list_filter[n].orbit_number] +=1
                        break
    
    return GDVMetric(GDV,vertex)

def graph_creation_from_XML(path):
    G = nx.Graph()   # or DiGraph, MultiGraph, MultiDiGraph, etc
    xmldoc = minidom.parse(path)
    itemlist = xmldoc.getElementsByTagName('edge') 
    for items in itemlist:
        src=items.attributes['source'].value
        trg=items.attributes['target'].value
        G.add_edge(int(src[1:]),int(trg[1:]))
    return G

def main():
    path = r'C:\Users\lohit\Desktop\RA\Datasets\naive_reduce\naive_reduce\naive_reduce\msg_size_1\without_ninja\run_6\slices\slice_8.graphml'
    graph = graph_creation_from_XML(path)
    nodes_list_graph = list(graph.nodes())
    GDV_matrix = []
    for node in nodes_list_graph:
        gdv_out = GDV_calculation(node,graph)
        GDV_matrix.append(gdv_out.GDV)
    print(np.matrix(GDV_matrix))


    
if __name__ == '__main__':
    main()






