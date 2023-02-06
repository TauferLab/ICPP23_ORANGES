#include "../include/orbitmatch.hpp"

int main()
{
    OrbitMatcher orbitmatcher = OrbitMatcher();
    // read_orbit_table(orbitmatcher);
    Edge edges[] = {{2,3}, {2,4}, {3,5}, {3,6}, {4,5}, {5,6}};
    std::vector<Edge> v_edges = {{2,3}, {2,4}, {3,5}, {3,6}, {4,5}, {5,6}};
    std::vector<Node> nodes = {2,3,4,5,6};
    std::vector<Node> endnodes;
    Motif motif;
    motif.root = 2;
    motif.nodes = nodes;
    motif.endnodes = endnodes;
    motif.edges = v_edges;
    std::cout << "motif created" << std::endl;
    //const Motif c_motif (motif);
    //std::make_shared<Edge>({2,3})
    
    //std::map<Node,std::array<int, 5>> orbits;
    orbitmatcher.get_orbits(motif);
    //for (auto const& [node, gdv]: orbits)
}