"""
    Serial Python implementation of the Fido algorithm for subtree counting
    It is customizable via the global variables/flags
    This file also has a bunch of testing and analysis available too
    This file only generates subtrees
    Author: Ali Khan
"""
from collections import namedtuple, Counter
import random
import ast
from itertools import combinations, chain
from timeit import default_timer
import networkx as nx


####### FLAGS ########
MAX_NODES = 5 # max number of nodes in the graphlet
# MAX_ORBIT = 72 # the maximum orbit number, gives orbits for all graphlet of MAX_NODES or less

TEST_GRAPH = nx.Graph()
# TEST_GRAPH.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,3),(2,4),(3,4),(4,5)])  # semi-house
# TEST_GRAPH.add_edges_from([(0,1),(0,2),(1,2)])  # triangle graph
# TEST_GRAPH.add_edges_from([(0,3),(3,4),(3,2),(2,5),(5,1),(4,1),(4,5)])  # house-porch
# TEST_GRAPH.add_edges_from([(0,1),(1,0),(2,1),(1,2),(3,2),(2,3),(3,0),(0,3),(4,0),(0,4),(4,1),(1,4)])  # house
# TEST_GRAPH = nx.tutte_graph()
# TEST_GRAPH = nx.petersen_graph()
# TEST_GRAPH = nx.complete_graph(15)
# TEST_GRAPH = nx.grid_2d_graph(4,4)  
TEST_GRAPH = nx.karate_club_graph()
# TEST_GRAPH = nx.gnp_random_graph(100,.2)

GRAPH_PERMUTATION = "KCore"  # Reorder the graph. Options: 'Degree', 'Random', 'None', 'Ego', or 'KCore'
EGO_RADIUS = 1  # radius for ego network. used only when using `Ego` graph permutation
DEBUG_SHOW_ORBITS = False
DEBUG_SHOW_REPEAT_TREES = True and DEBUG_SHOW_ORBITS # shows every nonunique tree. Do only if total number of subtrees is low
###### END FLAGS ######


####### Main algo ######
def enumerate_subtrees(G: nx.Graph) -> dict[int : nx.Graph]:
    """
    Calculates the GDV of each node
      Parameters:
        G (Networkx Graph): the graph to calculate the GDVs on
      Returns:
        dictionary of node: list of trees rooted at node
  """
    final_subtrees = {i:[] for i in G.nodes()}  # dict of root:list of subtrees found in graph
    # final_subtrees = {i:0 for i in G.nodes()}  # count of subtrees at each node
    subtree_queue = []

    time_vector = {}
    for root in G.nodes():
        start_time = default_timer()
        print("root:", root)

        # First build L(G) at level 1 with end nodes as neighbor of root
        neighbors = [u for u in G.neighbors(root) if u > root]
        # leaf_combinations = [combo for i in range(1, len(neighbors) + 1) for combo in combinations(neighbors, i)]
        leaf_combinations = [combo for i in range(1, MAX_NODES) for combo in combinations(neighbors, i)]
        for combo in leaf_combinations:
            if len(combo) + 1 <= MAX_NODES:
                T = nx.Graph()
                T.add_nodes_from(list((w for w in combo if w > root)),endnode=True)
                T.add_edges_from(list((root, w) for w in combo if w > root))
                subtree_queue.append(T)  # append the new extended subtree to L
                final_subtrees[root].append(T)
                # final_subtrees[root] += 1

        # Next, go through subtree queue and expand the subtree
        while subtree_queue:
            T = subtree_queue.pop(-1)
            if T.number_of_nodes() <= MAX_NODES:
                # end_nodes = [u for u, degree in T.degree() if degree == 1 and u != root]  # leaves of T
                end_nodes = [u for u,isendnode in T.nodes(data='endnode') if isendnode == True]
                visited = set(T.nodes())  # all nodes of T
                
                # get neighboring edges of T's endnodes in G
                neighboring_edges = []
                for v in end_nodes:
                    neighboring_edges.extend([(w,v) for w in G.neighbors(v) if w not in visited and w > root])
                
                # generate all combinations of the edges
                leaf_combinations = [combo for i in range(1, MAX_NODES - T.number_of_nodes() + 1) for combo in combinations(neighboring_edges, i)]
                # leaf_combinations = [combo for i in range(1, len(neighboring_edges) + 1) for combo in combinations(neighboring_edges, i)]

                for combo in leaf_combinations:
                    # generate new subtrees
                    if len(combo) + T.number_of_nodes() <= MAX_NODES:
                            # filter the combos. we dont want the combo if it uses the same node eg {(1,2),(1,3)}
                            edge_set = set([edge[0] for edge in combo])
                            if len(edge_set) != len(combo):
                                continue

                            # create new subtree
                            Tnew = nx.Graph()
                            Tnew.add_edges_from(T.edges)  # edges from its parent

                            # mark new endnodes
                            for edge in combo:
                                Tnew.add_node(edge[0],endnode=True)

                            Tnew.add_edges_from(combo)  # add the new edges to new subtree

                            subtree_queue.append(Tnew)  # append the new extended subtree to L
                            final_subtrees[root].append(Tnew)
                            # final_subtrees[root] += 1

        time_vector[root] = default_timer() - start_time  # time elapsed to process subtrees rooted at node i
    return final_subtrees, time_vector

###### END MAIN ALGO ######


###### RELABEL GRAPH ######
def _ego_relabeling(G, radius=1) -> list[int]:
    # relabel the graph based on the complexity of its ego graph
    ego_complexity = {}  
    for nodeid in G.nodes():
        _ego = nx.ego_graph(G,nodeid,radius=radius,undirected=True)
        ego_complexity[nodeid] = _ego.number_of_edges()
    print("Ego Complexity", ego_complexity)
    sorted_nodelist = list([x[0] for x in sorted(ego_complexity.items(),key=lambda kv: kv[1],reverse=False)])
    return sorted_nodelist    
    

def permute_graph(G, ordering, **kwargs)->nx.Graph:
    # relabel the graph. ordering can be None, Random, Degree, and Ego
    # 'Ego' also accepts the 'radius' of the ego network
    if ordering == "None":
        return G
    elif ordering == "Degree":
        sorted_nodelist = list([x[0] for x in sorted(G.degree(),key=lambda kv: kv[1],reverse=False)])  # lowest to highest degree
    elif ordering == "Ego":
        radius = kwargs.get("radius", 1)
        sorted_nodelist = _ego_relabeling(G, radius=radius)
    elif ordering == "KCore":
        sorted_nodelist = list([x[0] for x in sorted(nx.core_number(G).items(),key=lambda kv: kv[1])])
    elif ordering == "Random":
        sorted_nodelist = sorted(G.nodes(), key=lambda k: random.random())
    else:
        raise ValueError(f"Invalid argument for ordering: `{ordering}`")
    
    print("sorted_nodelist",sorted_nodelist)

    # create a mapping old label -> new label
    # node_mapping = dict(zip(G.nodes(), sorted(G.nodes(), key=lambda k: random.random()))) # random labeling
    node_mapping = {val:i for i,val in enumerate(sorted_nodelist)}
    print("Node Mapping", node_mapping)

    # build a new graph
    G = nx.relabel_nodes(G, node_mapping)
    return G

###### END RELABEL GRAPH ######

###### UTILS #######
def check_subtree_repeated(trees: list[nx.Graph]) -> list[tuple[tuple[int,int]]]:
    # takes list of subtrees and returns the subtrees that are repeated
    tree_set = set()
    repeats = []
    for tree in trees:
        sorted_tree = _normalize_edgelist(tree.edges())
        if sorted_tree in tree_set:
            repeats.append(sorted_tree)
        else:
            tree_set.add(sorted_tree)
    return repeats

def _normalize_edgelist(edges: list[int,int]) -> tuple[tuple[int,int]]:
    # normalize edgelist by sorting the edges
    sorted_edgelist = []
    # first sort the nodes in each edge
    for edge in edges:
        sorted_edge = tuple(edge)
        if edge[0] > edge[1]:
            sorted_edge = tuple([edge[1],edge[0]])
        sorted_edgelist.append(sorted_edge)
    # next, sort the edges
    sorted_edgelist.sort()
    return tuple(sorted_edgelist)

def flatten_tree(G:nx.Graph, node: int, seen: set) -> list[int]:
    # flatten tree into list. Preorder traversal
    # use: print(flatten_tree(T,root,set()))
    seen.add(node)
    out = [node]
    children = []
    for v in G.neighbors(node):
        if v not in seen:
            children.extend(flatten_tree(G,v,seen))
    if children:
      out.append(children)
    return out
###### END UTILS ######


def _test_subtree_repeats():
    T1 = nx.Graph()
    T1.add_edges_from([(3,4)])
    T2 = nx.Graph()
    T2.add_edges_from([(3,4)])
    T3 = nx.Graph()
    T3.add_edges_from([(0,1),(0,2)])
    T4 = nx.Graph()
    T4.add_edges_from([(0,2),(0,1)])
    repeated_subtree = check_subtree_repeated([tuple(T1.edges()),tuple(T2.edges()),tuple(T3.edges()),tuple(T4.edges())])
    print(repeated_subtree)
    print(len(repeated_subtree), "Repeated Trees")
    print(len(set(repeated_subtree)),"Uniquely repeated trees")

##### Orbit Matching #####

def get_orbit(motif: nx.Graph,root: int) -> int:
    # matches the motif to the correct graphlet orbit

    # check if the size of the motif is correct
    motif_n = motif.number_of_nodes()
    assert motif_n <= MAX_NODES

    # calculate degree signature for the motif
    motif_deg_sig = sorted(list(deg for node, deg in motif.degree(motif.nodes())))

    # calculate the distance signature for the motif
    motif_dist_sig = [0] * (motif_n)
    for dist,freq in Counter(dict(nx.shortest_path_length(motif,source=root)).values()).items():
        motif_dist_sig[dist] = freq

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

def get_orbit_info():
    # initialize the deg signature and distance signatures for graphlet orbits
    # orbitinfo.tsv only has trees upto 5 nodes

    OrbitInfo = namedtuple('OrbitInfo',['graphlet','orbit','deg_sig','dist_sig'])
    global TreeOrbitInfo
    TreeOrbitInfo = []
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


if __name__ == "__main__":

    # Preprocessing: relabel the nodes
    start = default_timer()
    G = TEST_GRAPH
    G = permute_graph(TEST_GRAPH,ordering=GRAPH_PERMUTATION, radius=EGO_RADIUS)
    print("Elapsed time for preprocessing:", default_timer() - start)

    # compute subtrees
    start = default_timer()
    subtrees_dict,times = enumerate_subtrees(G)
    # subtree_counts, times = enumerate_subtrees(G)
    endtime = default_timer()-start
    print("Time per node", dict(sorted(times.items(),key=lambda x: x[1])))

    # Count number of subtrees at node v in graph
    subtree_counts = {root:len(trees) for root,trees in sorted(subtrees_dict.items(),key=lambda x: x[0])}
    print("Subtree counts:",  subtree_counts)
    print("Number of subtrees", sum(subtree_counts.values()))
    print("Elapsed time for algo: ", endtime)

    # analyze the subtrees enumerate and count the repetitions
    subtrees = list(chain(*subtrees_dict.values()))
    repeated_subtree = check_subtree_repeated(subtrees)
    print(len(repeated_subtree), "Repeated Trees")
    print(len(set(repeated_subtree)),"Uniquely repeated trees")

    if len(repeated_subtree) != 0:
        get_orbit_info()  # initialize precalculated orbit info
        
        repeat_roots = {i:0 for i in G.nodes()}  # count of repeated subtrees rooted at node i
        normalized_subtrees = {root:[_normalize_edgelist(T.edges()) for T in tlist] for root,tlist in subtrees_dict.items()}
        if sum(subtree_counts.values()) < 100000:
            # takes too long for when number of subtrees is very large
            for i,tree in enumerate(repeated_subtree):
                # roots = []
                # find all the roots associated with the subtree
                for root,tree_list in normalized_subtrees.items():
                    # check if has multiple roots
                    # for tl in tree_list:
                    #     if tree == tl:
                    #         roots.append(root)
                    for tl in tree_list:
                        if tree == tl:
                            repeat_roots[root] += 1
                            if DEBUG_SHOW_ORBITS:
                                T = nx.Graph(incoming_graph_data=tree)
                                if DEBUG_SHOW_REPEAT_TREES:
                                    print(f"{i}\t Root {root} :: Orbit {get_orbit(T,root)} :: {flatten_tree(T, root, set())}")
                                else:
                                    print(f"{i})\t Root {root} :: Orbit {get_orbit(T,root)} :: {tree}")
                            break
            print("Repeats rooted at node v:", repeat_roots)

    # print("Printing all subtrees at root 1")
    # for i,tree in enumerate(subtrees_dict[1]):
    #     print(flatten_tree(tree, 1, set()))
