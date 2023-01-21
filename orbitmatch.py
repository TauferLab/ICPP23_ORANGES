import ast
from collections import namedtuple, Counter
import networkx as nx

MAX_NODES = 5 # max number of nodes in the graphlet

##### Orbit Matching #####
# compare the sorted degree vector and the shortestpath counts vector of the motif to the corresponding orbit

def get_orbit(motif) -> int:
    G = nx.Graph()
    G.add_edges_from(motif.edgelist())
    return get_orbit_nx(G, motif.root)

def get_orbit_nx(motif: nx.Graph, root: int) -> int:
    # matches the motif with the given root node to the correct graphlet orbit

    # check if the size of the motif is correct
    motif_n = motif.number_of_nodes()
    assert motif_n <= MAX_NODES

    # calculate degree signature for the motif
    motif_deg_sig = sorted(list(deg for node, deg in motif.degree(motif.nodes())))

    # calculate the distance signature for the motif
    motif_dist_sig = [0] * (motif_n)
    # for dist,freq in Counter(dict(nx.shortest_path_length(motif,source=root)).values()).items():
    #     motif_dist_sig[dist] = freq
    for u in motif.nodes():
        _sp = list(nx.all_shortest_paths(motif,root,u))
        motif_dist_sig[len(_sp[0])-1] += len(_sp)

    # match the rooted motif to the orbit
    for orbit_info in TreeOrbitInfo:
        orbit_n = len(orbit_info.deg_sig)
        if orbit_n == motif.number_of_nodes() \
            and orbit_info.deg_sig == motif_deg_sig \
            and orbit_info.dist_sig == motif_dist_sig:

            # print(f"Subgraph with degsig:{motif_deg_sig} and {motif_dist_sig} matches to orbit {orbit_info.orbit}")
            return orbit_info.orbit
    else:
        raise Exception("Orbit for given motif not found")

def get_orbits(motif_edges) -> dict:
    # returns dict of root nodes and the associated orbit number
    motif = nx.Graph()
    motif.add_edges_from(motif_edges.edgelist())
    return get_orbits_nx(motif)

def get_orbits_nx(motif: nx.Graph) -> dict:
    # returns dict of root nodes and the associated orbit number
    
    # output dict, key=node in motif, value=orbit
    orbits = dict()

    # check if the size of the motif is correct
    motif_n = motif.number_of_nodes()
    assert motif_n <= MAX_NODES

    # calculate degree signature for the motif
    motif_deg_sig = sorted(list(deg for node, deg in motif.degree(motif.nodes())))

    # calculate the distance signature for the motif
    # shortest_paths = dict(nx.shortest_path_length(motif))

    for node in motif.nodes():
        # calculate degree signature for the motif
        motif_dist_sig = [0] * (motif_n)
        # for dist,freq in  Counter(shortest_paths[node].values()).items():
        #     motif_dist_sig[dist] = freq
        for u in motif.nodes():
            _sp = list(nx.all_shortest_paths(motif,node,u))
            motif_dist_sig[len(_sp[0])-1] += len(_sp)
        # print(node,motif_deg_sig, motif_dist_sig)
        # match the rooted motif to the orbit
        for orbit_info in TreeOrbitInfo:
            orbit_n = len(orbit_info.deg_sig)
            if orbit_n == motif_n \
                and orbit_info.deg_sig == motif_deg_sig \
                and orbit_info.dist_sig == motif_dist_sig:

                # print(f"Subgraph with degsig:{motif_deg_sig} and {motif_dist_sig} matches to orbit {orbit_info.orbit}")
                # return orbit_info.orbit
                orbits[node] = orbit_info.orbit

    if len(orbits) != motif_n:
        raise Exception(f"Orbit for given motif and root ({node}) not found")
    else:
        return orbits

def get_orbit_info():
    # initialize the deg signature and distance signatures for graphlet orbits
    # orbitinfo.tsv only has trees upto 5 nodes

    OrbitInfo = namedtuple('OrbitInfo',['graphlet','orbit','deg_sig','dist_sig'])
    global TreeOrbitInfo
    TreeOrbitInfo = []
    # with open("orbitinfo_trees.tsv") as file:
    with open("orbitinfo.tsv") as file:
        for line in file.readlines()[1:]:
            row = line.split()
            graphlet_id = int(row[0])
            orbit_id = int(row[1])
            deg_sig = sorted(ast.literal_eval(row[2]))
            dist_sig = ast.literal_eval(row[3])
            orbit_info = OrbitInfo(graphlet=graphlet_id, orbit=orbit_id, deg_sig=deg_sig, dist_sig=dist_sig)
            TreeOrbitInfo.append(orbit_info)

def _test_get_orbit():
    get_orbit_info()
    G = nx.Graph()
    G.add_edges_from([(0,1)])
    assert get_orbit(G,0) == 0
    G = nx.Graph()
    G.add_edges_from([(0,1),(1,2),(2,3),(2,4)])
    assert get_orbit(G,2) == 21
    G = nx.Graph()
    G.add_edges_from([(1, 2), (1, 5), (2, 3), (4, 5)])
    print(get_orbit(G,1))
###### END Orbit Matching #####