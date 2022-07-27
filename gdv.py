"""
    Serial Python implementation of the Fido algorithm for gdv calculation
    Author: Ali Khan
"""
import networkx as nx
from itertools import combinations
from collections import defaultdict, Counter
from ast import literal_eval

"""
Fido Algorithm. Developed by Ali, Ujwal, Lohit, and Dr. Bhowmick
Input : Graph G
Output : L(G)[List of trees that grow over time]
       M(G)[List of all trees with back edges]
For each vertex in the Graph G:
    Fix the vertex as root
    Build L(G),trees at level 1 with end nodes as neighbour of root
    For each tree in L(G):
        For each end node:
            Find neighbours that are not visited in the tree.
            Make all combinations of the edges and add each combination to tree
            Update list of end nodes
                End nodes are nodes with degree one.
            Mark all trees that can be considered as sub trees
                If the end node is less than the root, then tree is not considered as sub tree
            Add trees to L(G)
            Find whether root tree can be considered for back edges and add it to M(G) and calculate GDV .
                      Back edge is added between the lowest node in the cycle and its lowest neighbour.
        End For
        If the Tree is considered as sub tree, then calculate its GDV.
    End For
End For
"""

MAX_NODES = 4
MAX_ORBIT = 14  # the maximum orbit number, gives orbits for all graphlet of MAX_NODES or less


def calc_gdv_graph(G: nx.Graph) -> list[list[int]]:
    """
    Calculates the GDV of each node
      Parameters:
        G (Networkx Graph): the graph to calculate the GDVs on
      Returns:
        List of GDV (list of ints) of each node
  """
    L = []  # list of trees
    final_subtrees = set()  # list of indices of L that are to be considered as a subtree

    gdvs = defaultdict(lambda: [0]*(MAX_ORBIT + 1))  # graphlet degree vector for all nodes
    for root in G.nodes():
        print("root:", root)
        visited = {root}  # mark root as visited
        # build L(G) at level 1 with end nodes as neighbor of root
        _count = 0
        for v in G.neighbors(root):
            T = nx.Graph(incoming_graph_data=[(root, v)])  # single edge tree with one end as root and other as neighbor of root
            L.append(T)
            if v >= root:  # if the end node is less than the root, then tree is not considered as a subtree
                _count += 1
                final_subtrees.add(T)
        print(f"Added {G.number_of_nodes()} subtrees to L")
        print(f"Added {_count} number of trees to final_subtrees")

        i = 0  # index of L
        while i < len(L) - 1:
            i += 1
            T = L[i]
            print(f"L{i}={T}")
            if T.number_of_nodes() <= MAX_NODES:
                end_nodes = [u for u, d in T.degree() if d == 1]
                while end_nodes:
                    v = end_nodes.pop()
                    visited.add(v)  # mark node v as visited in T
                    print("Examining end node", v)
                    # neighbors of end_node in G that havent been visited yet
                    neighbors = [w for w in G.neighbors(v) if w not in visited]  # TODO: make sure that all parents are marked visited
                    # enumerate the set of combinations of neighbors connected to the end_node
                    leaf_combinations = [combo for i in range(1, len(neighbors) + 1) for combo in combinations(neighbors, i)]
                    # create new subtrees of each permutation
                    for combo in leaf_combinations:
                        if len(combo) + T.number_of_nodes() <= MAX_NODES:
                            Tnew = T.copy()  # new subtree
                            Tnew.add_edges_from(list((v, w) for w in combo))  # add node(s) to the tree
                            L.append(Tnew)  # append the new extended subtree to L
                            if all(w >= root for w in combo):
                                final_subtrees.add(Tnew)  # mark as final subtree if all end nodes are greater than root
                    # check if tree can be considered for back edges
                    tree_edges = set(T.edges())
                    tree_nodes = set(T.nodes())
                    for u in T.nodes():
                        # for each node in the tree get back edges in T
                        back_edges = []
                        for v in G.neighbors(u):
                            if (u, v) not in tree_edges:  # check if (u,v) is a backedge
                                if is_valid_backedge(T, u, v, tree_nodes):
                                    back_edges.append((u, v))
                        # generate all motifs (we dont need to save them)
                        back_edge_motifs = [combo for i in range(len(back_edges) + 1) for combo in
                                            combinations(back_edges, i)]
                        for combo in back_edge_motifs:
                            motif = T.copy()
                            motif.add_edges_from(list((v, w)) for w in combo)
                            orbit = get_orbit(motif)
                            gdvs[root][orbit] += 1

            # if tree is considered as sub tree, calculate gdv
            if T in final_subtrees:
                orbit = get_orbit(T)
                gdvs[root][orbit] += 1

    return dict(gdvs)


def is_valid_backedge(G: nx.Graph, u: int, v: int, exclude=None) -> bool:
    """
        Checks whether edge (u,v) is a valid backedge.
        Validity is determined by if (u,v) is the lowest node in the cycle and its lowest neighbor
        Arguments:
            G (nx.Graph): the graph to get the backedges from
            u (int): a node in the edge
            v (int): the other node in the edge
            exclude (List[int]): list of nodes to exclude
        Returns:
            Boolean whether (u,v) is a valid backedge
    """
    # use iterative DFS to get the cycle at node u, if it exists
    # NOTE: may want to replace with something more performant in the future
    cycle = []
    stack = [u]
    if not exclude:
        seen = set()
    else:
        seen = set(exclude)
    while stack:
        u = stack.pop()
        if u == v:
            break
        for v in G.nodes:
            if v not in seen:
                seen.add(v)
                stack.append(v)
                cycle.append(v)
    # TODO: use naive cycle finding algorithm
    # may replace later

    # check if u or v is the lowest node in the cycle
    # if one of them is the lowest, check if the other is the lowest neighbor
    if cycle and ((min(cycle) == u and min(G.neighbors) == v) \
            or (min(cycle) == v and min(G.neighbors) == u)):
        return True
    return False


def get_orbit(motif: nx.Graph) -> int:
    """
        Finds the orbit that matches the given motif by comparing the degree and distance 
            signatures of the graphlet associated with an orbit against that of the motif.
        Degree signature of a graph is the sorted list of degrees of each node in the graph
        Distance signature of a graph is the sorted list of the number of nodes at distance x from the orbit node
        Arguments:
            motif (nx.Graph): a graph where |V| <= N_MOTIF_MAX
        Returns:
            The orbit number (int). Raises exception if no orbit matches
    """
    n = motif.number_of_nodes()
    assert n <= MAX_NODES
    # the degree signature of the motif, sorted degrees for
    degree_sig = sorted(list(deg for node, deg in motif.degree(motif.nodes())))
    # the distance signatures the list of sorted list of shortest paths rooted at i

    # dist_sigs = [list(Counter(list(distances.values())).values()) for u, distances in nx.shortest_path_length(motif)]
    dist_sigs = []  # the distance signatures of the motif, a list of the number of shortest paths x distance away
    for _, distances in nx.shortest_path_length(motif):
        counts = Counter(distances.values())
        dist_sig = list(counts.values())
        dist_sig.extend([0]*(6 - len(dist_sig)))  # TODO dont hardcode 6 (fix orbit_list.txt)
        dist_sigs.append(dist_sig)

    # iterate over possible orbits and return the one that has the matching deg sig and dist sig
    for orbit, signatures in orbit_sigs.items():
        orbit_deg_sig = sorted(signatures[0])
        orbit_dist_sig = signatures[1]
        num_nodes_graphlet = len(orbit_deg_sig)
        # compare the number of nodes and the degree and dist signatures
        if n == num_nodes_graphlet \
                and degree_sig == orbit_deg_sig \
                and orbit_dist_sig in dist_sigs:
            ## print(f"Subgraph {motif} maps to orbit {orbit}")
            return orbit
    else:  # if the loop ends, no orbit matched
        raise Exception("Orbit for given motif not found")


def gen_orbit_sig() -> dict[int, tuple[list[int], list[int]]]:
    """
        Generate and store the orbits, degree signature, and distance signatures associated 
        Returns:
            orbit_sigs, a dictionary of orbits and a tuple of the associated degree signature and distance signature for the graphlet
                where degree_signature (List[int]): the (sorted) degree signature for motif associated with the given orbit
                and distance_signature (List[int]): the (sorted) shortest path distances from the orbit node to every other
    """
    # TODO: get it to work for up to orbit 72 (all graphlets of 5 nodes)
    # TODO: Replace this code with algorithm that actually generates them instead of just reading hardcoded
    orbit_sigs = {}
    with open("./GRAAL/orbit_list.txt", 'r') as orbit_file:
        for line in orbit_file.readlines():
            orbit, deg_sig, dist_sig = line.split('/')
            orbit = int(orbit)
            deg_sig = literal_eval(deg_sig)
            dist_sig = literal_eval(dist_sig)
            orbit_sigs[orbit] = deg_sig, dist_sig
    return orbit_sigs


orbit_sigs = gen_orbit_sig()
