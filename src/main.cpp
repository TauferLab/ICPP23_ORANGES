/**
 * @file main.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-09-13
 * 
 * @copyright Copyright (c) 2022
 * 
 */
//INPUT HEADERS
#include "translate_from_input.hpp"
#include "input_to_network.hpp"
#include"structure_defs.hpp"

//OUTPUT HEADERS
#include "printout_network.hpp"
#include "printout_others.hpp"

// // HEADERS
// #include "CMPB_structure.hpp"
// #include  "subgraph_extract.hpp"
// #include "createN_Struct.hpp"

#include "ADJ/degree_centrality.hpp"

// #include "add_multiple_edge.hpp"
#include "ADJ/create_network.hpp"
/*** All Headers Required From ESSENS **/

#include <ctime>
#include <math.h>
#include <omp.h>
#include <algorithm>

int main(int argc, char *argv[])
{
    /***** Preprocessing the Graph ***********/
    clock_t q = clock();
    //Check if valid input is given
    if ( argc < 3) 
    { 
        std::cout << "INPUT ERROR:: At least 2 inputs required. First: filename \n Second: Filetypes: 1:node_node_wt 2:node_wt_node 3:node_node 4:node_node (Only option 1 is active now) \n Third. Name of new file \n Fourth. Name of Map file" << std::endl;; 
        return 0;
    }
    //Check to see if file opening succeeded
    std::ifstream the_file ( argv[1] ); if (!the_file.is_open() ) { cout<<"INPUT ERROR:: Could not open file\n";}
    
    A_Network X;
    int nodes=-1;
    map_int_st revmap;
    
    //Proces File if Option Given by
    if(argc==5)
    {
        // Preprocess Nodes to Numbers
        //Stores file in argv[3]: store map in argv[4]
        //Vertices start from 0
        int type=atoi(argv[2]);
        translate_input(argv[1],type,argv[3],argv[4]);
        
        //Remove Duplicate Edges and Self Loops; Create Undirected Graphs
        // process_to_simple_undirected();
        q = clock() - q;
        cout << "Total Time for Preprocessing"<< ((float)q)/CLOCKS_PER_SEC <<"\n";

        /******* Read Graph (GUI) and Create Reverse Map*****************/
        //Obtain the list of edges.
        q = clock();
        readin_network(&X,argv[3],nodes);
        
        //Create Reversemap
        
        nodes = X.size();
        create_map(argv[4],&revmap);
        
        q = clock() - q;
        cout << "Total Time for Reading Network" << ((float)q)/CLOCKS_PER_SEC <<"\n";
        /**** Read Graph (GUI) ***********/
    } else {
        
    }

    std::cout << "Finished" << std::endl;
}