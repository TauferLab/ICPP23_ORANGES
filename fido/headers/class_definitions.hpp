#ifndef CLASS_DEFINITIONS_HPP
#define CLASS_DEFINITIONS_HPP

#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

using namespace std;

class OrbitMetric
{    
  public: 
    vector<int> orbitDegree;
    vector<int> orbitDistance;
    int orbitNumber;
    OrbitMetric(int cOrbitNumber, vector<int> cOrbitDegree, vector<int> cOrbitDistance)
    {
      orbitNumber= cOrbitNumber;
      orbitDistance = cOrbitDistance;
      orbitDegree = cOrbitDegree;

    }
};

class GDVMetric{
public: 
  vector<int> GDV;
  int node;
  GDVMetric() {};
  GDVMetric(int cNode, vector<int> cGDV)
  {
    GDV = cGDV;
    node = cNode;
  };
  GDVMetric(const GDVMetric &old_metric) {
    node = old_metric.node;
    GDV = old_metric.GDV;
  }
};

class Orbits {
  public:
    Kokkos::View<int**> degree;
    Kokkos::View<int**> distance;
    Kokkos::View<int*> start_indices;

  Orbits() {}

  Orbits(int num_orbits, int max_nodes) {
    degree = Kokkos::View<int**>("Degree signatures", num_orbits, max_nodes);
    distance = Kokkos::View<int**>("Distance signatures", num_orbits, max_nodes+1);
    start_indices = Kokkos::View<int*>("Starting indices for graphlets", max_nodes+2);
  }

  int num_orbits() {
    return degree.extent(0);
  }
};

using Ordinal       = int;
using Offset        = int;
using Scalar        = float;
using device_type   = typename Kokkos::Device<Kokkos::DefaultExecutionSpace, 
                      typename Kokkos::DefaultExecutionSpace::memory_space>;
using matrix_type   = typename KokkosSparse::CrsMatrix<Scalar, Ordinal, device_type, void, Offset>;
using graph_type    = typename matrix_type::staticcrsgraph_type;
using row_map_type  = typename graph_type::row_map_type;
using entries_type  = typename graph_type::entries_type;
using values_type   = typename matrix_type::values_type;

using GDVs = Kokkos::View<int**>;

#endif
