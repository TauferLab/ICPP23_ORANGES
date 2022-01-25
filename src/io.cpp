#include "io.hpp"

//This method converts a string containing integers to a vector of integers
void convert_string_vector_int(string* str, vector<int>* output,string delimiter)
{
  size_t pos = 0;
  int token;
  string s;
  s = *str;
  do
  {
    pos = s.find(delimiter);
    token = stoi(s.substr(0, pos));
    output->push_back(token);
    s.erase(0, pos + delimiter.length());
  }
  while (pos != std::string::npos);

}

//This method takes the file and converts it into orbits and saves in output
void kokkos_readin_orbits(ifstream *file, Orbits& orbits )
{
  string line;
  string signature_delimiter;
  string internal_delimiter;
  signature_delimiter = "/";
  internal_delimiter= ",";
  int num_orbits = 0;
  file->clear();
  file->seekg(0);
  while(std::getline(*file,line))
  {
    num_orbits++;
  }
  file->clear();
  file->seekg(0);
  orbits = Orbits(num_orbits, 5);
  orbits.start_indices_host(0) = 0;
  orbits.start_indices_host(1) = 0;
  orbits.start_indices_host(2) = 0;
  size_t orbit_size = 2;
  int orbit_counter = 0;
  while(std::getline(*file,line))
  {
    string s= line;
    size_t pos = 0;
    vector<vector<int>> vector_line;
    do
    {
      vector<int> segment; 
      string token;
      pos = s.find(signature_delimiter);
      token = s.substr(0, pos);
      token.erase(remove(token.begin(), token.end(), '['), token.end());
      token.erase(remove(token.begin(), token.end(), ']'), token.end());
      convert_string_vector_int(&token,&segment,internal_delimiter);
      s.erase(0, pos + signature_delimiter.length());
      vector_line.push_back(segment);
    }
    while (pos!= std::string::npos);
    sort(vector_line[1].begin(), vector_line[1].end());
    if(vector_line[1].size() > orbit_size) {
      orbit_size = vector_line[1].size();
      orbits.start_indices_host(orbit_size) = orbit_counter;
    }
    for(size_t i=0; i<vector_line[1].size(); i++) {
      orbits.degree_host(orbit_counter, i) = vector_line[1][i];
    }
    for(size_t i=0; i<vector_line[2].size(); i++) {
      orbits.distance_host(orbit_counter, i) = vector_line[2][i];
    }
    orbit_counter++;
  }
  orbits.start_indices_host(6) = num_orbits;
  Kokkos::deep_copy(orbits.degree, orbits.degree_host);
  Kokkos::deep_copy(orbits.distance, orbits.distance_host);
  Kokkos::deep_copy(orbits.start_indices, orbits.start_indices_host);
}

bool sort_by_leading_edge(vector<int>& a, vector<int>& b) {
  return a[0] < b[0];
}

void readin_graph(ifstream* file, matrix_type& graph) 
{
  string line;;

  int lastrow = -1;
  vector<int> rows, cols;
  vector<float> vals;
  vector<vector<int>> edge_list;
  while(std::getline(*file,line))
  {
    int u, v;
    float w;
    sscanf(line.c_str(), "%d %d %f", &u, &v, &w);
    if(u > lastrow) 
      lastrow = u;
    if(v > lastrow) 
      lastrow = v;
    vector<int> edge1;
    edge1.push_back(u);
    edge1.push_back(v);
    if(w == 0.0 || w == 0) {
        edge1.push_back(0);
    } else {
        edge1.push_back(1);
    }
    edge_list.push_back(edge1);
    if(u != v) {
      vector<int> edge2;
      edge2.push_back(v);
      edge2.push_back(u);
      if(w == 0.0 || w == 0) {
            edge2.push_back(0);
      } else {
            edge2.push_back(1);
      }
      edge_list.push_back(edge2);
    }
  }
  sort(edge_list.begin(), edge_list.end(), sort_by_leading_edge);

  vector<int> rowmap(lastrow+2, 0);
  for(size_t i=0; i<edge_list.size(); i++) {
    rows.push_back(edge_list[i][0]);
    cols.push_back(edge_list[i][1]);
    vals.push_back(edge_list[i][2]);
  }
  
  for(size_t i=0; i<rows.size(); i++) {
    if(rows[i] != cols[i])
    {
      rowmap[rows[i]+1]++;
    }
    else
    {
      rowmap[rows[i]+1]++;
    }
  }
  for(size_t i=1; i<rowmap.size(); i++) {
    rowmap[i] += rowmap[i-1];
  }
  for(size_t i=0; i<rowmap.size()-1; i++) {
    sort(cols.begin()+rowmap[i], cols.begin()+rowmap[i+1]);
  }

  Kokkos::View<float*> vals_view("Vals", vals.size());
  Kokkos::View<int*> cols_view("Cols", cols.size());
  Kokkos::View<int*> row_map_view("Rows", rowmap.size());
  Kokkos::View<float*>::HostMirror vals_host = Kokkos::create_mirror_view(vals_view);
  Kokkos::View<int*>::HostMirror cols_host = Kokkos::create_mirror_view(cols_view);
  Kokkos::View<int*>::HostMirror row_map_host = Kokkos::create_mirror_view(row_map_view);
  for(int i=0; i<vals.size(); i++) 
    vals_host(i) = vals[i];
  for(int i=0; i<cols.size(); i++) 
    cols_host(i) = cols[i];
  for(int i=0; i<rowmap.size(); i++) 
    row_map_host(i) = rowmap[i];
  Kokkos::deep_copy(vals_view, vals_host);
  Kokkos::deep_copy(cols_view, cols_host);
  Kokkos::deep_copy(row_map_view, row_map_host);
  
  graph = matrix_type("Graph", lastrow+1, lastrow+1, vals.size(), 
                      vals_view, row_map_view, cols_view);
  return;
}

void write_similarity_matrix(Kokkos::View<double**>& sim_mat) {
  ofstream myfile; 
  Kokkos::View<double**>::HostMirror sim_mat_host = Kokkos::create_mirror_view(sim_mat);
  Kokkos::deep_copy(sim_mat_host, sim_mat);
  string filename = "out_similarity_matrix.txt";
  myfile.open(filename);
  for(int i=1; i<sim_mat.extent(0);i++) {
    for(int j=1;j<sim_mat.extent(1);j++) {
      myfile<<" { "<<sim_mat_host(i,j)<<" } ";
    }
    myfile<<"||"<<endl;
  }
  myfile.close();
}

void write_gdvs(GDVs::HostMirror& graph1_GDV_host, std::string& graph_tag1, 
                GDVs::HostMirror& graph2_GDV_host, std::string& graph_tag2) {
  // Print out GDVs into files
  ofstream myfile; 
  string gdv_file1 = "out_gdv_1_" + graph_tag1 + ".txt";
  string gdv_file2 = "out_gdv_2_" + graph_tag2 + ".txt";
  myfile.open(gdv_file1, ofstream::trunc);
  for (int i = 0; i < graph1_GDV_host.extent(0); i++) {
    for (int j = 0; j< graph1_GDV_host.extent(1); j++) {
      myfile << graph1_GDV_host(i,j) << " ";
    }
    myfile << endl;
  }
  myfile.close();
  myfile.open(gdv_file2, ofstream::trunc);
  for (int i = 0; i < graph2_GDV_host.extent(0); i++) {
    for (int j = 0; j< graph2_GDV_host.extent(1); j++) {
      myfile << graph2_GDV_host(i,j) << " ";
    }
    myfile << endl;
  }
  myfile.close();
}

