#ifndef PRINT_DISCONNECTED_GRAPH
#define PRINT_DISCONNECTED_GRAPH

#include "structure_defs.hpp"
#include <iostream>

void print_disconnected_network(A_Network &network) 
{
  for (int i = 0; i < network.size(); i++) {
    cout << "-------------------------------------" << endl;
    cout << network[i].Row << ":" << endl;
    for (int j = 0; j < network[i].ListW.size(); j++) {
      cout << network[i].ListW[j].first << " , " << endl;
    }
    cout << "-------------------------------------" << endl;
  }
  return;
}

#endif
