#from networkit import *
#import networkit as nk
import igraph
import os
import argparse
import re


def main ( graph_txt, graphml_file_path ):

    G = igraph.Graph.Read_Ncol(graph_txt, names=True)
    #G = nk.readGraph(graph_txt, nk.Format.GraphML)
    if not os.path.isdir(graphml_file_path):
        os.makedirs(graphml_file_path)
    #H = nk.Graph(G)
    filename = graphml_file_path + graph_txt.split('/')[len(graph_txt.split('/'))-1].split('.')[0] + '.graphml'
    print(filename)
    G.write_graphml(filename)
    #nk.writeGraph(H, filename, nk.Format.METIS)

if __name__ == "__main__":
    desc = "Converts a comm pattern graphml file to a metis file."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("graph_txt",
                        help="Path to comm pattern file")
    parser.add_argument("graphml_file_path",
                        help="Path to output where metis file will be stored")
    args = parser.parse_args()
    main( args.graph_txt, args.graphml_file_path )
