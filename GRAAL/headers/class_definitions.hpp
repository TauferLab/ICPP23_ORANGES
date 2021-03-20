

#ifndef CLASS_DEFINITIONS_HPP
#define CLASS_DEFINITIONS_HPP

using namespace std;

class OrbitMetric
{
    public:
    vector<int> orbitDegree;
    vector<int> orbitDistance;
    int orbitNumber;
    OrbitMetric(int cOrbitNumber, vector<int> cOrbitDegree, vector<int> cOrbitDistance)
    {
        orbitNumber   = cOrbitNumber;
        orbitDistance = cOrbitDistance;
        orbitDegree   = cOrbitDegree;
    }
};

class GDVMetric
{
    public:
    vector<int> GDV;
    int node;
    GDVMetric() {};
    GDVMetric(int cNode, vector<int> cGDV)
    {
        GDV  = cGDV;
        node = cNode;
    };
    GDVMetric(const GDVMetric& old_metric)
    {
        node = old_metric.node;
        GDV  = old_metric.GDV;
    }
};

#endif
