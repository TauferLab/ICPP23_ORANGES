// Required directives are to be added from ESSENS to use the functionality

#include "ADJ/find_Xneighbors.hpp"
#include "RAJA/RAJA.hpp"
#include "RAJA/policy/cuda/policy.hpp"
#include "RAJA/util/sort.hpp"
#include <umpire/ResourceManager.hpp>
#include "headers/GDV_functions.hpp"
#include "headers/GPUGDV_functions.hpp"
#include "headers/class_definitions.hpp"
#include "headers/print_disconnected_graph.hpp"
#include "input_to_network.hpp"
#include "printout_network.hpp"
#include "printout_others.hpp"
#include "structure_defs.hpp"
#include <math.h>
#include <time.h>

void Calculate_GDV(int, A_Network, vector<OrbitMetric>&, GDVMetric&);
void readin_orbits(ifstream*, vector<OrbitMetric>*);
void convert_string_vector_int(string*, vector<int>*, string);
void GDV_vector_calculation(A_Network graph, vector<GDVMetric>* graph_GDV,
        vector<OrbitMetric> orbits, int p);
void metric_formula(GDVMetric gdvm, double* gdv_score);
double GDV_distance_calculation(GDVMetric gdvm1, GDVMetric gdvm2);
void Similarity_Metric_calculation_for_two_graphs(A_Network graph1, A_Network graph2,
        vector<OrbitMetric> orbits, int p);
using namespace std;
using inner_policy = RAJA::cuda_exec<256>;
using outer_policy = RAJA::seq_exec;

int main(int argc, char* argv[])
{
    clock_t tStart = clock();
    clock_t q, q1, q2, t;
    GDV_functions gdvf;
    /* Accepts file as input.
       Input should be in the format of ( node1 node2 weight ).
       By default, weight should be 1
       Should be sorted. */
    ifstream the_file1(argv[1]);
    if (!the_file1.is_open())
    {
        cout << "INPUT ERROR:: Could not open the graph input file\n";
    }

    ifstream the_file2(argv[2]);
    if (!the_file2.is_open())
    {
        cout << "INPUT ERROR:: Could not open the second graph input file\n";
    }

    ifstream the_file3(argv[3]);
    if (!the_file3.is_open())
    {
        cout << "INPUT ERROR:: Could not open the orbit file\n";
    }
    int p = atoi(argv[4]);
    vector<OrbitMetric> orbits;
    readin_orbits(&the_file3, &orbits);
    // Objects for testing orbit creation
    // print_vector(orbits[1].orbitDegree);
    // vector<OrbitMetric> filter_o = gdvf.orbit_filter(&orbits,3);
    /* Read in Network: Reads the file and converts it into a network of type A_Network*/
    A_Network X;
    readin_network(&X, argv[1], 0, -1);


    A_Network Y;
    readin_network(&Y, argv[2], 0, -1);
    GDV_functions test_gdvf;

    // objects for testing orbit list
    // vector<OrbitMetric> filtered_orbits ;
    // gdvf.orbit_filter(orbits,3,filtered_orbits);
    // print_vector(filtered_orbits[1].orbitDistance);

    // Objects for testing GDV induced subgraph function
    // A_Network subgraph;
    // vector<int> subgraph_nodes;
    // subgraph_nodes.push_back(0);
    // subgraph_nodes.push_back(1);
    // subgraph_nodes.push_back(2);
    // subgraph_nodes.push_back(7);
    // gdvf.inducedSubgraph(X, subgraph_nodes, subgraph);
    // print_disconnected_network(subgraph);
    // cout<<"subgraph for 0,1,2,6"<<endl;
    // print_network(subgraph);

    // Objects for testing connectedness function
    // bool is_connected = false;
    // gdvf.isConnected(subgraph, is_connected);
    // cout << is_connected << endl;

    // Objects for testing degree signature
    // vector<int> degree_sig;
    // test_gdvf.degree_signature(X, degree_sig);
    // print_vector(degree_sig);

    // Objects for testing distance signature
    // vector<int> distance_sig;
    // test_gdvf.distance_signature(2, X, distance_sig);
    // print_vector(distance_sig);
    // for (int i:X)
    // {
    //   // Calculate_GDV(i,X);
    // }

    Similarity_Metric_calculation_for_two_graphs(X, Y, orbits, p);
    printf("Time taken: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);
    return 0;
}

void Similarity_Metric_calculation_for_two_graphs(A_Network graph1, A_Network graph2,
        vector<OrbitMetric> orbits, int p)
{
    vector<GDVMetric> graph1_GDV;
    vector<GDVMetric> graph2_GDV;
    // double start = omp_get_wtime();
    GDV_vector_calculation(graph1, &graph1_GDV, orbits, p);
    GDV_vector_calculation(graph2, &graph2_GDV, orbits, p);

    // vector<vector<int>> similarity_matrix( graph1_GDV.size() , vector<int> (graph2_GDV.size(),
    // 0));
    int m = (int) graph1_GDV.size();
    int n = (int) graph2_GDV.size();

    double sim_mat[m][n];
    for (GDVMetric gdvm1 : graph1_GDV)
    {
        for (GDVMetric gdvm2 : graph2_GDV)
        {
            sim_mat[gdvm1.node][gdvm2.node] = GDV_distance_calculation(gdvm1, gdvm2);
        }
    }
    // double end = omp_get_wtime();
    // cout<<"Omp time taken is " <<end - start<<endl;
    ofstream myfile;
    string filename;
    filename = "out_similarity_matrix.txt";
    myfile.open(filename);

    for (int i = 1; i < m; i++)
    {
        for (int j = 1; j < n; j++)
        {
            myfile << " { " << sim_mat[i][j] << " } ";
        }
        myfile << "||" << endl;
    }
    myfile.close();
}

double GDV_distance_calculation(GDVMetric gdvm1, GDVMetric gdvm2)
{
    double gdv1_score;
    double gdv2_score;
    metric_formula(gdvm1, &gdv1_score);
    metric_formula(gdvm2, &gdv2_score);

    double similarity_score = abs(gdv1_score - gdv2_score);

    return similarity_score;
}

void metric_formula(GDVMetric gdvm, double* gdv_score)
{
    int sum = 0;
    //  Formula used here is the vector norm 2 or the l2 norm.
    //  We square each value and do summation. Now take a square root.
    for (int x : gdvm.GDV)
    {
        int number;
        number = x * x;
        sum    = sum + number;
    }
    *gdv_score = sqrt(sum);
}

void GDV_vector_calculation(A_Network graph, vector<GDVMetric>* graph_GDV,
        vector<OrbitMetric> orbits, int p)
{
    // omp_set_num_threads(p);
    //#pragma omp parallel for num_threads(p) schedule(static)
    RAJA::forall<outer_policy>(RAJA::RangeSegment(0, graph.size()), [&](int i) {
            ADJ_Bundle node = graph[i];
            vector<int> GDV_1;
            GDVMetric gdvMetric(node.Row, GDV_1);
            Calculate_GDV(node.Row, graph, orbits, gdvMetric);
            // cout<<"gdv for node "<<node.Row<<endl;
            graph_GDV->push_back(gdvMetric);
            });
    // for (ADJ_Bundle node:graph)
    // {
    // vector<int> GDV_1;
    // GDVMetric gdvMetric(node.Row,GDV_1);
    // Calculate_GDV(node.Row,graph,orbits,gdvMetric);
    // // cout<<"gdv for node "<<node.Row<<endl;
    // graph_GDV->push_back(gdvMetric);
    // }
}
void Calculate_GDV(int node, A_Network Graph, vector<OrbitMetric>& orbits, GDVMetric& gdvMetric)
{
    auto& res = umpire::ResourceManager::getInstance();
    umpire::Allocator halloc = res.getAllocator("HOST");
    umpire::Allocator dalloc = res.getAllocator("DEVICE");
    GDV_functions gdvf;
    GPUGDV_functions ggdvf;
    intvecvec combinationsList;
    A_Network_raw raw_Graph;
    orbvec raw_orbits;

    stdnet_to_rawnet(Graph, raw_Graph, halloc);
    stdorb_to_raworb(orbits, raw_orbits, halloc);
    //vector<int> gdv(orbits.size(),0);
    intvec gdv = new_intvec_umpire(orbits.size(), dalloc);
    res.memset(gdv.vec, 0, orbits.size());
    // printf("calculating GDV for node %d\n",node);
    vector<int> neighbours;

    gdvf.find_neighbours(node, Graph, 4, &neighbours);
    // print_vector(neighbours);
    int set[neighbours.size()];
    std::copy(neighbours.begin(), neighbours.end(), set);
    int numElements       = *(&set + 1) - set;

    int combinations_size = 0;
    for (int node_count = 1; node_count < 5; node_count++){
        combinations_size += ((ggdvf.fact(neighbours.size()))/(ggdvf.fact(node_count)*ggdvf.fact(neighbours.size()-node_count)));
    }
    combinationsList = new_intvecvec_umpire(combinations_size, halloc);
    for (int node_count = 1; node_count < 5; node_count++)
    {
        ggdvf.find_combinations_raw(&set[0], numElements, node_count, combinationsList, halloc);
    }

    for(int i = 0; i < combinations_size; i++){
        combinationsList.vec[i].vec = (int*) res.move(combinationsList.vec[i].vec, dalloc);
    }
    combinationsList.vec = (intvec*) res.move(combinationsList.vec, dalloc);

    for(int i = 0; i < raw_Graph.nodes_len; i++){
        raw_Graph.vec[i].ListW.vec = (Edge_raw*) res.move(raw_Graph.vec[i].ListW.vec, dalloc);
        //raw_Graph.vec[i].Ops.vec = (int*) res.move(raw_Graph.vec[i].Ops.vec, dalloc);
    }
    raw_Graph.vec = (Adjlist*) res.move(raw_Graph.vec, dalloc);

    for(int i = 0; i < raw_orbits.veclen; i++){
        raw_orbits.vec[i].orbitDegree.vec = (int*) res.move(raw_orbits.vec[i].orbitDegree.vec, dalloc);       
        raw_orbits.vec[i].orbitDistance.vec = (int*) res.move(raw_orbits.vec[i].orbitDistance.vec, dalloc);
    }
    raw_orbits.vec = (OrbitMetric_raw*) res.move(raw_orbits.vec, dalloc);

    // cout<<"Node count is "<<node_count<<endl;
    // cout<<"total combinations are : "<<combinationsList.size()<<endl;
    //for (vector<int> combination : combinationsList)

    RAJA::forall<inner_policy>(RAJA::RangeSegment(0,combinations_size), [node, raw_orbits, raw_Graph, combinationsList, gdv, ggdvf] FIDO_DEVICE (int i){
            A_Network_raw induced_sgraph = new_network(6);

            //vector<int> subgraph_degree_signature;
            //vector<int> subgraph_distance_signature;
            intvec subgraph_degree_signature = new_intvec(6);
            intvec subgraph_distance_signature = new_intvec(6);
            bool is_connected = false;
            //combination.push_back(node);
            pushback_intvec(combinationsList.vec[i], node);
            ggdvf.inducedSubgraph_raw(raw_Graph, combinationsList.vec[i], induced_sgraph);
            ggdvf.isConnected_raw(induced_sgraph, is_connected);
            if (is_connected)
            {
            ggdvf.degree_signature_raw(induced_sgraph, subgraph_degree_signature);
            ggdvf.distance_signature_raw(node, induced_sgraph, subgraph_distance_signature);
            //vector<OrbitMetric> filter_orbits;
            orbvec filter_orbits = new_orbvec(22);
            ggdvf.orbit_filter_raw(raw_orbits, combinationsList.vec[i].veclen, filter_orbits);
            //for (OrbitMetric orbit : filter_orbits)
            for(int i = 0; i < filter_orbits.veclen; i++) 
            {
                OrbitMetric_raw orbit = filter_orbits.vec[i];
                //sort(orbit.orbitDegree.begin(), orbit.orbitDegree.end());
                //RAJA::sort<inner_policy>(orbit.orbitDegree.vec, orbit.orbitDegree.vec+orbit.orbitDegree.veclen);
                RAJA::detail::intro_sort(orbit.orbitDegree.vec, orbit.orbitDegree.vec+orbit.orbitDegree.veclen, 
                        [] FIDO_DEVICE (int a, int b){return a < b;}, 10);
                //sort(subgraph_degree_signature.begin(), subgraph_degree_signature.end());
                //RAJA::sort<inner_policy>(subgraph_degree_signature.vec, subgraph_degree_signature.vec+subgraph_degree_signature.veclen);
                RAJA::detail::intro_sort(subgraph_degree_signature.vec, subgraph_degree_signature.vec+subgraph_degree_signature.veclen,
                        [] FIDO_DEVICE (int a, int b){return a < b;}, 10);
                if (eq_intvec(orbit.orbitDistance, subgraph_distance_signature) &&
                        eq_intvec(orbit.orbitDegree, subgraph_degree_signature))
                {
                    gdv.vec[orbit.orbitNumber] += 1;
                    break;
                }
            }
            delete_orbvec(filter_orbits);
            }
            delete_intvec(subgraph_degree_signature);
            delete_intvec(subgraph_distance_signature);
            for(int i = 0; i < induced_sgraph.nodes_len; i++){
                delete_adjlist(induced_sgraph.vec[i]);
            }
            delete_network(induced_sgraph);
    });
    gdv.vec = (int*) res.move(gdv.vec, halloc);
    std::copy(gdv.vec, gdv.vec+gdv.veclen, gdvMetric.GDV.begin());
    //gdvMetric.GDV  = gdv;
    gdvMetric.node = node;
    delete_intvec_umpire(gdv, dalloc);
    for(int i = 0; i < combinations_size; i++){
        delete_intvec_umpire(combinationsList.vec[i], dalloc);
    }
    delete_intvecvec_umpire(combinationsList, dalloc);
    for(int i = 0; i < raw_Graph.nodes_len; i++){
        delete_adjlist_umpire(raw_Graph.vec[i], dalloc);
    }
    delete_network_umpire(raw_Graph, dalloc);

    for(int i = 0; i < raw_orbits.veclen; i++){
        delete_orbitmetric_umpire(raw_orbits.vec[i], dalloc);
    }
    delete_orbvec_umpire(raw_orbits, dalloc);
}

// This method takes the file and converts it into orbits and saves in output
void readin_orbits(ifstream* file, vector<OrbitMetric>* output)
{
    string line;
    string signature_delimiter;
    string internal_delimiter;
    signature_delimiter = "/";
    internal_delimiter  = ",";
    while (std::getline(*file, line))
    {
        string s   = line;
        size_t pos = 0;
        vector<vector<int>> vector_line;
        do
        {
            vector<int> segment;
            string token;
            pos   = s.find(signature_delimiter);
            token = s.substr(0, pos);
            token.erase(remove(token.begin(), token.end(), '['), token.end());
            token.erase(remove(token.begin(), token.end(), ']'), token.end());
            convert_string_vector_int(&token, &segment, internal_delimiter);
            s.erase(0, pos + signature_delimiter.length());
            vector_line.push_back(segment);
        } while (pos != std::string::npos);

        OrbitMetric orbMetric(vector_line[0][0], vector_line[1], vector_line[2]);
        output->push_back(orbMetric);
    }
}

// This method converts a string containing integers to a vector of integers
void convert_string_vector_int(string* str, vector<int>* output, string delimiter)
{
    size_t pos = 0;
    int token;
    string s;
    s = *str;
    do
    {
        pos   = s.find(delimiter);
        token = stoi(s.substr(0, pos));
        output->push_back(token);
        s.erase(0, pos + delimiter.length());
    } while (pos != std::string::npos);
}
