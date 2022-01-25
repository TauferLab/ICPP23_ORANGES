#ifndef IO_HPP
#define IO_HPP

#include <fstream>
#include <string>
#include <vector>
#include "class_definitions.hpp"
#include <Kokkos_Core.hpp>

//This method takes the file and converts it into orbits and saves in output
void kokkos_readin_orbits(ifstream *file, Orbits& orbits );

void readin_graph(ifstream* file, matrix_type& graph);

void write_similarity_matrix(Kokkos::View<double**>& sim_mat);

void write_gdvs(GDVs::HostMirror& graph1_GDV, std::string& graph_tag1, GDVs::HostMirror& graph2_GDV, std::string& graph_tab2);

#endif
