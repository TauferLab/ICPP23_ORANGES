#ifndef KOKKOS_FUNCTIONS_HPP
#define KOKKOS_FUNCTIONS_HPP

#include <vector>
#include <iterator>
#include <iostream>
#include <stdio.h>
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

using namespace std;
KOKKOS_INLINE_FUNCTION int64_t factorial(int n) {
  int64_t result = 1;
  for(int i=1; i<n; i++) {
    result *= i;
  }
  return result;
}

KOKKOS_INLINE_FUNCTION uint64_t get_num_combinations(int n, int k) {
  uint64_t result = 1;
  for(int i=n-k+1; i<=n; i++) {
    result *= i;
  }
  for(int i=2; i<=k; i++) {
    result /= i;
  }
  return result;
}

KOKKOS_INLINE_FUNCTION uint32_t
get_approx_num_neighbors(const matrix_type& graph, int node, int distance) {
  uint32_t num_neighbors = 0;
  auto row = graph.row(node);
  if(distance > 1) {
    for(size_t i=0; i<row.length; i++) {
      num_neighbors += get_approx_num_neighbors(graph, row.colidx(i), distance-1);
    }
  } else {
    return row.length;
  }
  return num_neighbors;
}

template<class IndexView>
KOKKOS_INLINE_FUNCTION void combination_from_position(IndexView& indices, int64_t position, int size, int k) {
  int64_t m = position;
  int n = size;
  int r = k;
  for(int i=0; i<indices.extent(0); i++) {
    indices(i) = 0;
  }
//  Kokkos::deep_copy(indices, 0);
  int64_t y;
  while(n > 0) {
    if(n > r && r >= 0) {
      y = get_num_combinations(n-1, r);
    } else {
      y = 0;
    }
    if(m >= y) {
      m -= y;
      indices(n-1) = 1;
      r -= 1;
    } else {
      indices(n-1) = 0;
    }
    n -= 1;
  }
  for(size_t i=0; i<k; i++) {
    int index = 0;
    for(size_t j=i; j<indices.extent(0); j++) {
      if(indices(j) == 1) {
        index = j;
        indices(j) = 0;
        break;
      }
    }
    indices(i) = index;
  }
  return;
}

void transpose_gdvs(GDVs& src, Kokkos::View<uint32_t**, Kokkos::LayoutLeft>& dst) {
  Kokkos::MDRangePolicy<Kokkos::Rank<2>> mdrange({0,0}, {src.extent(0), src.extent(1)});
  Kokkos::parallel_for("Transpose", mdrange, KOKKOS_LAMBDA(const int x, const int y) {
    dst(y,x) = src(x,y);
  });
}

template<class IndexView, class SetView, class CombinationView>
KOKKOS_INLINE_FUNCTION void kokkos_get_combination(IndexView& indices, 
                                                    const SetView& set, 
                                                    int k,
                                                    CombinationView& combo)
{
  for(int i=0; i<k; i++) {
    if(i < set.extent(0))
      combo(i) = set(indices(i));
  }
  return;
}

//KOKKOS_INLINE_FUNCTION Kokkos::View<uint32_t**, Kokkos::LayoutLeft> make_layout_left(GDVs gdv) {
//  Kokkos::View<uint32_t**, Kokkos::LayoutLeft> new_gdv("new gdv", gdv.extent(0), gdv.extent(1));
//  if(std::is_same<Kokkos::LayoutLeft,  GDVs::array_layout>::value) {
////    return gdv;
//  } else {
//    for(int i=0; i<gdv.extent(0); i++) {
//      for(int j=0; j<gdv.extent(1); j++) {
//        new_gdv(i,j) = gdv(i,j);
//      }
//    }
//    return new_gdv;
//  }
//}


namespace EssensKokkos {

  template <class SubviewType>
  KOKKOS_INLINE_FUNCTION int
  get_num_neighbors(const matrix_type& graph, int node, int distance, SubviewType& neighbors) {
    int num_neighbors = 0;
    int level = 0;
    int level_start = 0;
    int level_end = 0;
    auto root_row = graph.row(node);
    for(int i=0; i<root_row.length; i++) {
      if(root_row.colidx(i) > node) {
        neighbors(num_neighbors) = root_row.colidx(i);
        num_neighbors++;
      }
    }
    level_end = num_neighbors;
    level = 1;
    while(level < distance) {
      for(int i=level_start; i<level_end; i++) {
        auto row = graph.row(neighbors(i));
        for(int j=0; j<row.length; j++) {
          bool new_node = true;
          int u = row.colidx(j);
          if(u > node) {
            for(int k=0; k<num_neighbors; k++) {
              if(u == neighbors(k)) {
                new_node = false;
                break;
              }
            }
            if(new_node) {
              neighbors(num_neighbors) = u;
              num_neighbors++;
            }
          }
        }
      }
      level++;
      level_start = level_end;
      level_end = num_neighbors;
    }
    return num_neighbors;
  }

  template <class SubviewType>
  KOKKOS_INLINE_FUNCTION int
//  find_neighbours(int node, const matrix_type& graph,int distance, Kokkos::View<int*>& neighbors)
  find_neighbours(int node, const matrix_type& graph,int distance, SubviewType& neighbors)
  {
    int num_neighbors = 0;
    int level = 0;
    int level_start = 0;
    int level_end = 0;
    auto root_row = graph.row(node);
    for(int i=0; i<root_row.length; i++) {
      if(root_row.colidx(i) > node) {
        neighbors(num_neighbors) = root_row.colidx(i);
        num_neighbors++;
      }
    }
    level_end = num_neighbors;
    level = 1;
    while(level < distance) {
      for(int i=level_start; i<level_end; i++) {
        auto row = graph.row(neighbors(i));
        for(int j=0; j<row.length; j++) {
          bool new_node = true;
          int u = row.colidx(j);
          if(u > node) {
            for(int k=0; k<num_neighbors; k++) {
              if(u == neighbors(k)) {
                new_node = false;
                break;
              }
            }
            if(new_node) {
              neighbors(num_neighbors) = u;
              num_neighbors++;
            }
          }
        }
      }
      level++;
      level_start = level_end;
      level_end = num_neighbors;
    }
    return num_neighbors;
  }

  KOKKOS_INLINE_FUNCTION int
  get_num_subgraph_edges(const matrix_type& graph, const Kokkos::View<int*>& nodes)
  {
    int edge_count = 0;
    for(int i=0; i<nodes.size(); i++) {
      auto row = graph.rowConst(nodes(i));
      for(int j=0; j<row.length; j++) {
        for(int k=0; k<nodes.size(); k++) {
          if(row.colidx(j) == nodes(k)) {
            edge_count++;
          }
        }
      }
    }
    return edge_count;
  }

//  // Vector nodes should correspond to actual node labels, not node indices.
//  KOKKOS_INLINE_FUNCTION void 
////  kokkos_inducedSubgraph(const matrix_type& graph, const Kokkos::View<int*> &nodes, Kokkos::View<int*>& rowmap, Kokkos::View<int*> cols, Kokkos::View<float*> vals, matrix_type& output)
//  kokkos_inducedSubgraph(const matrix_type& graph, const Kokkos::View<int*> &nodes, matrix_type& output)
//  {
////cout << "Nodes: ";
////for(int i=0; i<nodes.size(); i++) {
////  cout << nodes(i) << " ";
////}
////cout << endl;
//    int edge_count = 0;
//    for(int i=0; i<nodes.size(); i++) {
//      auto row = graph.rowConst(nodes(i));
//      for(int j=0; j<row.length; j++) {
//        for(int k=0; k<nodes.size(); k++) {
//          if(row.colidx(j) == nodes(k))
//            edge_count++;
//        }
//      }
//    }
////cout << "Number of edges in subgraph: " << edge_count << endl;
//    Kokkos::View<int*> rowmap("Induced subgraph rowmap", nodes.size()+1);
//    Kokkos::View<int*> cols("Subgraph columns", edge_count);
//    Kokkos::View<float*> vals("Subgraph values", edge_count);
//    int edge_counter = 0;
//    for(int i=0; i<nodes.size(); i++) {
//      auto row = graph.row(nodes(i));
//      for(int j=0; j<row.length; j++) {
//        for(int k=0; k<nodes.size(); k++) {
//          if(row.colidx(j) == nodes(k) && (k != i)) {
//            cols(edge_counter) = k;
//            vals(edge_counter) = row.value(j);
//            rowmap(i+1)++;
//            edge_counter++;
//          }
//        }
//      }
//    }
//    rowmap(0) = 0;
//    for(int i=1; i<nodes.size()+1; i++) {
//      rowmap(i) += rowmap(i-1);
//    }
//    output = matrix_type("Induced subgraph", nodes.size(), nodes.size(), vals.size(), 
//                          vals, rowmap, cols);
////                          vals.data(), rowmap.data(), cols.data());
//    return;
//  }

  // Vector nodes should correspond to actual node labels, not node indices.
  template<class RowmapView, class CombinationView, class ColumnView, class ValueView>
  KOKKOS_INLINE_FUNCTION int
  kokkos_inducedSubgraph(const matrix_type& graph, const CombinationView& nodes, RowmapView& rowmap, ColumnView& cols, ValueView& vals)
  {
    int edge_count = 0;
    for(int i=0; i<nodes.size(); i++) {
      auto row = graph.rowConst(nodes(i));
      for(int j=0; j<row.length; j++) {
        for(int k=0; k<nodes.size(); k++) {
          if(row.colidx(j) == nodes(k))
            edge_count++;
        }
      }
    }
    int edge_counter = 0;
    for(int i=0; i<nodes.size(); i++) {
      auto row = graph.row(nodes(i));
      for(int j=0; j<row.length; j++) {
        for(int k=0; k<nodes.size(); k++) {
          if(row.colidx(j) == nodes(k) && (k != i)) {
            cols(edge_counter) = k;
            vals(edge_counter) = row.value(j);
            rowmap(i+1)++;
            edge_counter++;
          }
        }
      }
    }
    rowmap(0) = 0;
    for(int i=1; i<nodes.size()+1; i++) {
      rowmap(i) += rowmap(i-1);
    }
    return edge_count;;
  }

  // Vector nodes should correspond to actual node labels, not node indices.
  template<class CombinationView, class SubgraphType>
  KOKKOS_INLINE_FUNCTION int
  kokkos_induced_subgraph(const matrix_type& graph, const CombinationView& nodes, SubgraphType& subgraph)
  {
    for(size_t i=0; i<subgraph.extent(0); i++) {
      for(size_t j=0; j<subgraph.extent(1); j++) {
        subgraph(i,j) = 0;
      }
    }
    int edge_count = 0;
    for(size_t i=0; i<nodes.size(); i++) {
      auto row = graph.row(nodes(i));
      for(size_t j=0; j<row.length; j++) {
        for(size_t k=0; k<nodes.size(); k++) {
          if(row.colidx(j) == nodes(k) && (k!=i)) {
            subgraph(i,k) = row.value(j);
            edge_count++;
          }
        }
      }
    }
    return edge_count;
  }

//  // Check if a graph is connected with BFS.
//  KOKKOS_INLINE_FUNCTION
//  bool isConnected(const matrix_type& graph)
//  {
//    int connected_nodes = 0;
//    Kokkos::View<bool*> visited("Visited nodes", graph.numRows());
//    Kokkos::View<int*> queue("BFS queue", graph.numRows());
//    int queue_length = 0;
//    visited(0) = true;
//    queue(0) = 0;
//    queue_length = 1;
//    connected_nodes++;
//    while(queue_length > 0) {
//      int node = queue(queue_length-1);
//      queue_length -= 1;
//      auto row = graph.row(node);
//      for(int i=0; i<row.length; i++) {
//        if(!visited(row.colidx(i))) {
//          visited(row.colidx(i)) = true;
//          queue(queue_length) = row.colidx(i);
//          queue_length++;
//          connected_nodes++;
//        }
//      }
//    }
//    if(connected_nodes == graph.numRows()) {
//      return true;
//    }
//    return false;
//  }
  
//  // Check if a graph is connected with BFS.
//  template<class GraphType, class VisitedView, class QueueView>
//  KOKKOS_INLINE_FUNCTION
//  bool isConnected(const GraphType& graph, VisitedView& visited, QueueView& queue)
//  {
//    for(int i=0; i<visited.size(); i++) {
//      visited(i) = false;
//      queue(i) = -1;
//    }
//    int connected_nodes = 0;
//    int queue_length = 0;
//    visited(0) = true;
//    queue(0) = 0;
//    queue_length = 1;
//    connected_nodes++;
//    while(queue_length > 0) {
//      int node = queue(queue_length-1);
//      queue_length -= 1;
//      auto row = graph.row(node);
//      for(int i=0; i<row.length; i++) {
//        if(!visited(row.colidx(i))) {
//          visited(row.colidx(i)) = true;
//          queue(queue_length) = row.colidx(i);
//          queue_length++;
//          connected_nodes++;
//        }
//      }
//    }
//    if(connected_nodes == graph.numRows()) {
//      return true;
//    }
//    return false;
//  }
  
  // Check if a graph is connected with BFS.
  template<class GraphType, class VisitedView, class QueueView>
  KOKKOS_INLINE_FUNCTION
  bool is_connected(const GraphType& graph, VisitedView& visited, QueueView& queue)
  {
    for(size_t i=0; i<visited.size(); i++) {
      visited(i) = false;
      queue(i) = -1;
    }
    size_t connected_nodes = 0;
    int queue_length = 0;
    visited(0) = true;
    queue(0) = 0;
    queue_length = 1;
    connected_nodes++;
    while(queue_length > 0) {
      int node = queue(queue_length-1);
      queue_length -= 1;
      for(size_t i=0; i<graph.extent(1); i++) {
        if(!visited(i) && graph(node,i) == 1) {
          visited(i) = true;
          queue(queue_length) = i;
          queue_length++;
          connected_nodes++;
        }
      }
    }
    if(connected_nodes == graph.extent(0)) {
      return true;
    }
    return false;
  }
  
  // Calculate degree signature for network
  template<class GraphType, class ViewType>
  KOKKOS_INLINE_FUNCTION
  void calc_degree_signature(const GraphType& graph, ViewType& deg_sig)
  {
    if(deg_sig.size() != graph.extent(0)) {
      printf("Warning: Degree signature View does not match # of subgraph nodes.\n");
    } else {
      for(size_t i=0; i<deg_sig.size(); i++) {
        deg_sig(i) = 0;
      }
      for(size_t i=0; i<graph.extent(0); i++) {
        for(size_t j=0; j<graph.extent(1); j++) {
          if(graph(i,j) == 1)
            deg_sig(i) += 1;
        }
      }
      // Sort
      for(size_t i=0; i<deg_sig.size()-1; i++) {
        for(size_t j=0; j<deg_sig.size()-1-i; j++) {
          if(deg_sig(j) > deg_sig(j+1)) {
            int temp = deg_sig(j);
            deg_sig(j) = deg_sig(j+1);
            deg_sig(j+1) = temp;
          }
        }
      }
    }
    return;
  }
  
  // Calculate degree signature for network
  template<class GraphType, class ViewType>
  KOKKOS_INLINE_FUNCTION void 
  degree_signature(const GraphType& graph, ViewType& deg_sig)
  {
    if(deg_sig.size() != graph.numRows()) {
      std::cout << "Warning: Degree signature View does not match # of subgraph nodes." << endl;
    } else {
      for(int i=0; i<graph.numRows(); i++) {
        deg_sig(i) = graph.row(i).length; 
      }
    }
    for(int i=0; i<deg_sig.size()-1; i++) {
      for(int j=0; j<deg_sig.size()-1-i; j++) {
        if(deg_sig(j) > deg_sig(j+1)) {
          int temp = deg_sig(j);
          deg_sig(j) = deg_sig(j+1);
          deg_sig(j+1) = temp;
        }
      }
    }
    return;
  }


//  // Calculate distance signature for network
//  template<class ViewType>
//  KOKKOS_INLINE_FUNCTION
//  void distance_signature(int node, const matrix_type& graph, ViewType& dist_sig)
//  {
//    Kokkos::View<int*> queue("BFS queue", graph.numRows());
//    Kokkos::View<int*> distance("Distance View", graph.numRows());
//    Kokkos::View<bool*> visited("Visited nodes", graph.numRows());
//    Kokkos::deep_copy(distance, INT_MAX);
////    Kokkos::deep_copy(dist_sig, 0);
//    for(int i=0; i<dist_sig.size(); i++) {
//      dist_sig(i) = 0;
//    }
//    int queue_length = 0;
//    visited(node) = true;
//    queue(0) = node;
//    queue_length = 1;
//    distance(node) = 0;
//    while(queue_length > 0) {
//      int u = queue(queue_length-1);
//      queue_length -= 1;
//      auto row = graph.row(u);
//      for(int i=0; i<row.length; i++) {
//        if(!visited(row.colidx(i))) {
//          visited(row.colidx(i)) = true;
//          if(distance(u)+row.value(i) < distance(row.colidx(i))) {
//            distance(row.colidx(i)) = distance(u) + row.value(i);
//          }
//          queue(queue_length) = row.colidx(i);
//          queue_length++;
//        }
//      }
//    }
////cout << "Distance View: " ;
////for(int i=0; i<distance.size(); i++) {
////  cout << distance(i) << " ";
////}
////cout << endl;
//    int numrows = graph.numRows();
//    for(int i=0; i<numrows; i++) {
//      dist_sig(distance(i)) += 1;
//    }
//    return;
//  }

  // Calculate distance signature for network
  template<class GraphType, class ViewType, class VisitedView, class QueueView, class DistanceView>
  KOKKOS_INLINE_FUNCTION
  void distance_signature(int node, const GraphType& graph, ViewType& dist_sig, VisitedView& visited, QueueView& queue, DistanceView& distance)
  {
    for(int i=0; i<distance.size(); i++) {
      distance(i) = INT_MAX;
    }
    for(int i=0; i<visited.size(); i++) {
      visited(i) = false;
    }
    for(int i=0; i<dist_sig.size(); i++) {
      dist_sig(i) = 0;
    }
    int queue_length = 0;
    visited(node) = true;
    queue(0) = node;
    queue_length = 1;
    distance(node) = 0;
    while(queue_length > 0) {
      int u = queue(queue_length-1);
      queue_length -= 1;
      auto row = graph.row(u);
      for(int i=0; i<row.length; i++) {
        if(!visited(row.colidx(i))) {
          visited(row.colidx(i)) = true;
          if(distance(u)+row.value(i) < distance(row.colidx(i))) {
            distance(row.colidx(i)) = distance(u) + row.value(i);
          }
          queue(queue_length) = row.colidx(i);
          queue_length++;
        }
      }
    }
    int numrows = graph.numRows();
    for(int i=0; i<numrows; i++) {
      dist_sig(distance(i)) += 1;
    }
    return;
  }

  // Calculate distance signature for network
  template<class GraphType, class ViewType, class VisitedView, class QueueView, class DistanceView>
  KOKKOS_INLINE_FUNCTION
  void calc_distance_signature(int node, const GraphType& graph, ViewType& dist_sig, VisitedView& visited, QueueView& queue, DistanceView& distance)
  {
    for(size_t i=0; i<distance.size(); i++) {
      distance(i) = INT_MAX;
    }
    for(size_t i=0; i<visited.size(); i++) {
      visited(i) = false;
    }
    for(size_t i=0; i<dist_sig.size(); i++) {
      dist_sig(i) = 0;
    }
    int numcols = graph.extent(1);
    int queue_length = 0;
    visited(node) = true;
    queue(0) = node;
    queue_length = 1;
    distance(node) = 0;
    while(queue_length > 0) {
      int u = queue(queue_length-1);
      queue_length -= 1;
      for(int i=0; i<numcols; i++) {
        if(visited(i) == false && graph(u,i) == 1) {
          queue(queue_length) = i;
          queue_length++;
          visited(i) = true;
        }
        if(graph(u,i) == 1) {
          if(distance(i) > distance(u)+graph(u,i)) {
            dist_sig(distance(u)+graph(u,i)) += 1;
            if(distance(i) != INT_MAX)
              dist_sig(distance(i)) -= 1;
            distance(i) = distance(u) + graph(u,i);
          } else if(distance(i) == distance(u) + graph(u,i)) {
            dist_sig(distance(i)) += 1;
          }
        }
//        if(!visited(i) && graph(u,i) == 1) {
//          visited(i) = true;
//          if(distance(u)+graph(u,i) < distance(i)) {
//            distance(i) = distance(u) + graph(u,i);
//          }
//          queue(queue_length) = i;
//          queue_length++;
//        }
      }
    }
    dist_sig(0) = 1;
//    for(int i=0; i<numrows; i++) {
//      dist_sig(distance(i)) += 1;
//    }
    return;
  }

  template<class SignatureA, class SignatureB>
  KOKKOS_INLINE_FUNCTION
  bool compare_signatures(const SignatureA& sig1, const SignatureB& sig2) {
    for(size_t i=0; i<sig1.size(); i++) {
      if(sig1(i) != sig2(i)) 
        return false;
    }
    return true;
  }
}
#endif
