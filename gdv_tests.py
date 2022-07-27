import gdv
import networkx as nx
import csv


def petersen():
    G = nx.petersen_graph()
    gdvs = gdv.calc_gdv_graph(G)
    print(gdvs)
    with open("fido-orbit-counts-petersen.out", 'w', newline="") as f:
        csv.writer(f, delimiter=" ").writerows(gdvs.values())


if __name__ == "__main__":
    petersen()
