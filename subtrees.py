"""
    Serial Python implementation of the Fido algorithm for subtree counting
    It is customizable via the global variables/flags
    This file also has a bunch of testing and analysis available too
    This file only generates subtrees
    Author: Ali Khan
"""
import random
from itertools import combinations, chain
from timeit import default_timer
import json
import networkx as nx
from orbitmatch import *

####### FLAGS ########
RANDOM_SEED = 42
MAX_NODES = 5 # max number of nodes in the graphlet
# MAX_ORBIT = 72 # the maximum orbit number, gives orbits for all graphlet of MAX_NODES or less
CALCULATE_GDV = True  # otherwise only does subtrees
COUNTS_ONLY = True  # have the algorithm only produce the subtree counts instead of returning all the subtrees found

TEST_GRAPH = nx.Graph()
# TEST_GRAPH = nx.read_graphml("event_graph.graphml")
# TEST_GRAPH = nx.read_graphml("eventGraph1391.grapml")  # AMG med
# TEST_GRAPH = nx.read_graphml("eventgraph3523.graphml").to_undirected()  # AMG large
# print("radius",nx.radius(TEST_GRAPH), "diameter", nx.diameter(TEST_GRAPH))
# TEST_GRAPH = nx.read_graphml("event_graphMSG.graphml").to_undirected()
# TEST_GRAPH.add_edges_from([(0,1),(0,2),(1,2),(1,3),(2,3),(2,4),(3,4),(4,5)])  # semi-house
# TEST_GRAPH.add_edges_from([(0,1),(0,2),(1,2)])  # triangle graph
# TEST_GRAPH.add_edges_from([(0,3),(3,4),(3,2),(2,5),(5,1),(4,1),(4,5)])  # house-porch
# TEST_GRAPH.add_edges_from([(0,1),(1,0),(2,1),(1,2),(3,2),(2,3),(3,0),(0,3),(4,0),(0,4),(4,1),(1,4)])  # house
# TEST_GRAPH = nx.complete_graph(10)
# TEST_GRAPH = nx.tutte_graph()
# TEST_GRAPH = nx.petersen_graph()
# TEST_GRAPH = nx.complete_graph(15)
# TEST_GRAPH = nx.grid_2d_graph(4,4)  
TEST_GRAPH = nx.karate_club_graph()
# TEST_GRAPH = nx.read_gml('dolphins.gml')
# import scipy as sp 
# import scipy.io
# TEST_GRAPH = nx.from_scipy_sparse_matrix(sp.io.mmread("soc-wiki-Vote.mtx"))
# TEST_GRAPH = nx.read_edgelist("Wiki-Vote.txt")
# TEST_GRAPH = nx.gnp_random_graph(100,.2, seed=RANDOM_SEED, directed=False)
# TEST_GRAPH = nx.barabasi_albert_graph(100, 99, seed=RANDOM_SEED)
GRAPH_ORDERING = "KCore"  # Reorder the graph. Options: 'Degree', 'Random', 'None', 'Ego', or 'KCore'
EGO_RADIUS = 1  # radius for ego network. used only when using `Ego` graph permutation
DEBUG_SHOW_ORBITS = not COUNTS_ONLY and False
DEBUG_SHOW_REPEAT_TREES = True and DEBUG_SHOW_ORBITS # shows every nonunique tree. Do only if total number of subtrees is low
###### END FLAGS ######


####### Main algo ######
def enumerate_subtrees(G: nx.Graph) -> dict[int : nx.Graph]:
    """
    Calculates the GDV of each node
      Parameters:
        G (Networkx Graph): the graph to calculate the GDVs on
      Returns:
        dictionary of node: either a dictionary of trees rooted at node or a dictionary of the counts of trees rooted at each node
  """
    if COUNTS_ONLY:
        final_subtrees = {i:0 for i in G.nodes()}  # count of subtrees at each node
    else:
        final_subtrees = {i:[] for i in sorted(G.nodes())}  # dict of root:list of subtrees found in graph

    subtree_queue = []
    tree_orbits = {0:0,1:1,2:1,4:3,5:3,6:4,7:4,15:9,16:9,17:9,18:10,19:10,20:10,21:10,22:11,23:11}  # maps orbit to corresponding tree graphlet
    tree_counts = {x:0 for x in set(tree_orbits.values())}  # tracks how many tree graphlets occur in the graph
    orbit_counts = {x:0 for x in tree_orbits.keys()}
    queue_size = {}  # track the maximum size of the queue

    time_vector = {}
    for root in G.nodes():
        start_time = default_timer()
        print("root:", root)
        
        # First build L(G) at level 1 with end nodes as neighbor of root
        neighbors = [u for u in G.neighbors(root) if u > root]
        leaf_combinations = [combo for i in range(1, MAX_NODES) for combo in combinations(neighbors, i)]
        for combo in leaf_combinations:
            
            edges = list((root, w) for w in combo if w > root)
            nodes = list((w for w in combo if w > root))
            T = nx.Graph()
            T.add_nodes_from(nodes, endnode=True)  # add neighbors and mark endnodes
            T.add_edges_from(edges)        # add the edges

            subtree_queue.append(T)  # append the new extended subtree to L
            if COUNTS_ONLY:
                final_subtrees[root] += 1
            else:
                final_subtrees[root].append(T)
            
            for node, orbit in get_orbits_nx(T).items():
                orbit_counts[orbit] += 1

        queue_size[root] = len(subtree_queue)

        # Next, go through subtree queue and expand the subtree
        while subtree_queue:
            if queue_size[root] < len(subtree_queue):
                queue_size[root] = len(subtree_queue)   
            T = subtree_queue.pop(-1)
            if T.number_of_nodes() <= MAX_NODES:
                end_nodes = [u for u,isendnode in T.nodes(data='endnode') if isendnode == True]
                visited = set(T.nodes())  # all nodes of T
                
                # get neighboring edges of T's endnodes in G
                neighboring_edges = []
                for v in end_nodes:
                    neighboring_edges.extend([(w,v) for w in G.neighbors(v) if w not in visited and w > root])
                
                # generate all combinations of the edges
                leaf_combinations = [combo for i in range(1, MAX_NODES - T.number_of_nodes() + 1) for combo in combinations(neighboring_edges, i)]

                for combo in leaf_combinations:
                    # generate new subtrees
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

                    if COUNTS_ONLY:
                        final_subtrees[root] += 1
                    else:
                        final_subtrees[root].append(Tnew)
                    
                    for node, orbit in get_orbits_nx(T).items():
                        orbit_counts[orbit] += 1


        time_vector[root] = default_timer() - start_time  # time elapsed to process subtrees rooted at node i
    return final_subtrees, time_vector, queue_size, orbit_counts

###### END MAIN ALGO ######

###### BACKEDGES ######
def _test_backedges():
    T = nx.Graph()
    T.add_edges_from([(0,4),(2,4),(3,4)])
    backedges = [(2,3)]
    print("Tree:",flatten_tree(T,0,set()))
    print("Backedges:", backedges)
    # mindict = get_mindict(T)
    motifs = add_backedges(T,backedges)
    print("Motifs found:")
    for motif in motifs:
        print(f"\t{motif.edges()}")
    assert len(motifs) == 1

def get_mindict(T: nx.Graph) -> dict[tuple[int,int]:tuple[int,int]]:
    # return a dictionary of tuple of smallest 2 nodes in shortest path between two nodes, keyed by all pairs of nodes in the tree
    # TODO: can reduce size by half
    mindict= {}
    shortestpaths = dict(nx.all_pairs_shortest_path(T))
    sortednodes = sorted([node for node in T.nodes()],reverse=True)
    for node_a in sortednodes:
        for node_b in sortednodes:
            shortestpath = shortestpaths[node_a][node_b]
            sorted_shortestpaths = sorted([node for node in shortestpath])
            if len(sorted_shortestpaths) >= 2:
                mindict[(node_a,node_b)] = (sorted_shortestpaths[0],sorted_shortestpaths[1])
    return mindict

def add_backedges(T: nx.Graph, backedges: list[tuple[int,int]]) -> list[nx.Graph]:
    # add all backedge combinations to the tree if the backedge connects the two lowest nodes in the cycle
    # min_dict is a precomputed dict of tuples that contains the two minimum nodes in the shortest path from node_a to node_b
    valid_backedges = []
    out = []
    min_dict = get_mindict(T)
    for backedge in backedges:
        min1, min2 = min_dict[backedge]
        if min1 in backedge and min2 in backedge and min1 != min2:
            valid_backedges.append(backedge)
    backedge_combos = [combo for i in range(1,len(valid_backedges)+1) for combo in combinations(valid_backedges, i)]
    for backedge_combo in backedge_combos:
        Tnew = nx.Graph()
        Tnew.add_edges_from(T.edges())
        Tnew.add_edges_from(backedge_combo)
        out.append(Tnew)
    return out
###### END Backedges ######

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
    

def relabel_graph(G, ordering, **kwargs)->nx.Graph:
    # relabel the graph. ordering can be None, Random, Degree, and Ego
    # 'Ego' also accepts the 'radius' of the ego network
    print("Ordering: ", ordering)
    if ordering == "None":
        return G
    elif ordering == "Degree":
        sorted_nodelist = list([x[0] for x in sorted(G.degree(),key=lambda kv: kv[1],reverse=False)])  # lowest to highest degree
    elif ordering == "Ego":
        radius = kwargs.get("radius", 1)
        sorted_nodelist = _ego_relabeling(G, radius=radius)
    elif ordering == "KCore":
        print("Core numbers", nx.core_number(G))
        sorted_nodelist = list([x[0] for x in sorted(nx.core_number(G).items(),key=lambda kv: kv[1])])
    elif ordering == "Random":
        sorted_nodelist = sorted(G.nodes(), key=lambda k: random.random())
    else:
        raise ValueError(f"Invalid argument for ordering: `{ordering}`")
    
    # print("sorted_nodelist",sorted_nodelist)

    # create a mapping old label -> new label
    # node_mapping = dict(zip(G.nodes(), sorted(G.nodes(), key=lambda k: random.random()))) # random labeling
    node_mapping = {val:i for i,val in enumerate(sorted_nodelist)}
    # print("Node Mapping", node_mapping)

    # build a new graph
    G_relabel = nx.relabel_nodes(G, node_mapping, copy=True)
    G = nx.Graph()
    G.add_nodes_from(sorted(G_relabel.nodes))
    G.add_edges_from(G_relabel.edges())

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

def write_dict_to_file(fname, data):
    with open(fname,"w") as f:
        f.write(json.dumps(data))
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




if __name__ == "__main__":

    get_orbit_info()  # initialize precalculated orbit info

    # Preprocessing: relabel the nodes
    G = TEST_GRAPH
    print(G)
    start = default_timer()
    G = relabel_graph(TEST_GRAPH,ordering=GRAPH_ORDERING, radius=EGO_RADIUS)
    print("Elapsed time for preprocessing:", default_timer() - start)

    # compute subtrees
    start = default_timer()
    results, times, queue_sizes, tree_counts = enumerate_subtrees(G)
    endtime = default_timer()-start
    
    print("Time per node", times)
    write_dict_to_file("times_default.out",times)
    
    if COUNTS_ONLY:
        subtree_counts = results
        write_dict_to_file("counts.out",subtree_counts)
        
        # Count number of subtrees at node v in graph
        print("Subtree counts:",  subtree_counts)
        print("Maximal subtree queue sizes:", queue_sizes)
        print("Number of subtrees", sum(subtree_counts.values()))
        max_queue_node = max(queue_sizes,key=queue_sizes.get)
        print("Max queue size:", queue_sizes[max_queue_node], "  node:",max_queue_node)
        print("Elapsed time for algo: ", endtime)
        print("Tree Counts:",tree_counts)
        print(TEST_GRAPH)
    else:
        subtrees_dict = results
        subtree_counts = {root:len(trees) for root,trees in sorted(subtrees_dict.items(),key=lambda x: x[0])}
        write_dict_to_file("counts.out",subtree_counts)

        # Count number of subtrees at node v in graph
        print("Subtree counts:",  subtree_counts)
        print("Number of subtrees", sum(subtree_counts.values()))
        print("Elapsed time for algo: ", endtime)
        print("Tree Counts:",tree_counts)
        print(TEST_GRAPH)

        # analyze the subtrees enumerate and count the repetitions
        subtrees = list(chain(*subtrees_dict.values()))
        repeated_subtree = check_subtree_repeated(subtrees)
        print(len(repeated_subtree), "Repeated Trees")
        print(len(set(repeated_subtree)),"Uniquely repeated trees")

        if len(repeated_subtree) != 0:
            # get_orbit_info()  # initialize precalculated orbit info
            
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

