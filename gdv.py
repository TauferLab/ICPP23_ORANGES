"""
    Serial Python implementation of the Fido algorithm for gdv calculation
    Author: Ali Khan
"""
import networkx as nx
from itertools import permutations

"""
Fido Algorithm. Developed by Ujwal, Lohit, and Dr. Bhowmick
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
def calc_gdv_graph(G: nx.Graph) -> list[int]:
  """
    Calculates the GDV of each node
      Parameters:
        G (Networkx Graph): the graph to calculate the GDVs on
      Returns:
        List of GDV of each node
  """
  L = []  # list of trees
  M = []  # list of trees with backedges
  _final_subtree = []  # list of indices of L that are to be considered as a subtree
  gdvs = []
  for u in G.nodes():
    root = u
    # build L(G) at level 1 with end nodes as neighbor of root
    for v in G.neighbors(root):
        T = nx.Graph((root,v))  # single edge tree with one end as root and other as neighbor of root
        nx.set_node_attributes(T, {root:True, v:False}, name='visited')  # mark root as visited and neighbor as not visited
        # nx.set_node_attributes(T, {root:False, v:True}, name='end_node')  # mark neighbor as end-node
        L.append(T)
        if v >= root:  # if the end node is less than the root, then tree is not considered as a subtree
            _final_subtree.append(len(L)-1)

    i = 0  # index of L
    while i < len(L):
        i += 1
        T = L[i]
        # for i,T in enumerate(L):
        # end_nodes = [u for u, d in T.degree() if d == 1]
        end_nodes = [u for u, d in T.degree() if d == 1]
        while end_nodes:
            v = end_nodes.pop()
            T.nodes[v]['visited'] = True
            not_visited = [u for u, data in T.nodes(data=True) if not data['visited']]
            neighbors = [w for w in G.neighbors(v) if not G.nodes(w)['visited']]
            # permutations of leaves
            leaf_permutations = [combo for i in range(len(neighbors) + 1) for combo in permutations(neighbors, i) ]
            # create subtrees of each permutation
            for combo in leaf_permutations:
                Tnew = T.copy()  # copy of parent
                Tnew.add_edges(list((v,w) for w in combo))  # add node(s) to the tree
                L.append(Tnew)  # append the subtree to L
                if all(w >= root for w in combo):
                    _final_subtree.append(i)  # mark subtree if all end nodes are greater than root
            # check if tree can be considered for back edgesg
            # TODO  does the tree have to be a valid subtree?
            # TODO do we check all generated subtrees?
            # TODO where is this stored?
            # TODO does this get returned?
            tree_edges = set(T.edges())
            tree_nodes = set(T.nodes())
            for u in T.nodes():
                # for each node in the tree get back edges in T
                back_edges = []
                for v in G.neighbors[u]:
                    if (u,v) not in tree_edges:  # check if (u,v) is a backedge
                        if is_valid_backedge(T, u, v, tree_nodes):
                            back_edges.append((u,v))
                # generate all motifs
                back_edge_motifs = [combo for i in range(len(back_edges) + 1) for combo in permutations(back_edges, i) ]
                for combo in back_edge_motifs:
                    motif = T.copy()
                    motif.add_edge(list((v,w)) for w in combo)
                    M.append(gdv(motif, u))
                
        # if tree is considered as sub tree, calculate gdv
        # TODO how are we calculating the GDV?
        M.append(T)
        gdv(T, root)
  return M

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
    # basically DFS to get the cycle
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

    # check if u or v is the lowest node in the cycle
    # if one of them is the lowest, check if the other is the lowest neighbor
    if (min(cycle) == u and min(G.neighbors) == v) \
        or (min(cycle) == v and min(G.neighbors) == u):
        return True
    return False


def gdv(G: nx.Graph, node: int) -> int:
    """
        Calculates the GDV of a single node? in the graph
    """
    raise NotImplementedError

def gen_orbit_sig(orbit: int) -> tuple(int, list[int], list[int]):
    """
        Calculates and returns the degree signature and distance signature for the given graphlet orbit
        Arguments:
            orbit (int): the orbit number
        Returns:
            graphlet_id (int): the motif number associated with the given orbit
            degree_signature (List[int]): the (sorted) degree signature for motif associated with the given orbit
            distance_signature (List[int]): the (sorted) shortest path distances from the orbit node to every other
    """
    raise NotImplementedError