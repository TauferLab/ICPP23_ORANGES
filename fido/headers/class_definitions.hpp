#ifndef CLASS_DEFINITIONS_HPP
#define CLASS_DEFINITIONS_HPP

#include <Kokkos_Core.hpp>

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

  Orbits() {}

  Orbits(int num_orbits, int max_nodes) {
    degree = Kokkos::View<int**>("Degree signatures", num_orbits, max_nodes);
    distance = Kokkos::View<int**>("Distance signatures", num_orbits, max_nodes+1);
  }
};

#endif
