#from networkit import *
import networkit as nk
import os
import argparse
import re


def main ( comm_pattern, metis_file_path ):

    G = nk.readGraph(comm_pattern, nk.Format.GraphML)
    if not os.path.isdir(metis_file_path):
        os.makedirs(metis_file_path)
    H = nk.Graph(G)
    filename = metis_file_path + comm_pattern.split('/')[len(comm_pattern.split('/'))-1].split('.')[0] + '.graph'
    print(filename)
    nk.writeGraph(H, filename, nk.Format.METIS)

if __name__ == "__main__":
    desc = "Converts a comm pattern graphml file to a metis file."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("comm_pattern",
                        help="Path to comm pattern file")
    parser.add_argument("metis_file_path",
                        help="Path to output where metis file will be stored")
    args = parser.parse_args()
    main( args.comm_pattern, args.metis_file_path )
