#ifndef RAW_VECS_HPP
#define RAW_VECS_HPP 

// int_double
struct Edge
{
  int first;
  double second;
};

//vector<int>
struct intvec
{
  int* vec;
  int veclen;
};

//vector<int_double>
struct edgevec
{
  Edge* vec;
  int veclen;
};

//ADJ_Bundle
struct Adjlist
{
  int Row;
  edgevec ListW;
  intvec Ops;
  int list_len;
  int ops_len;
};

//vector<vector<int>>
struct intvecvec
{
   intvec* vec;
   int veclen;
}

//A_Network
struct A_Network_raw
{
  Adjlist* vec;
  int nodes_len;
};

//OrbitMetric
struct OrbitMetric_raw
{
  intvec orbitDegree;
  intvec orbitDistance;
  int orbitNumber;
};

//vector<OrbitMetric>
struct orbvec
{
  OrbitMetric_raw *vec;
  int veclen;
};

void pushback_edgevec(edgevec &vec, Edge in){
  vec->vec[vec->veclen++] = in;
}

void pushback_intvec(intvec &vec, int in){
  vec->vec[vec->veclen++] = in;
}

void pushback_adjlist(A_Network_raw &vec, Adjlist in){
  vec->vec[vec->veclen++] = in;
}

void pushback_intvecvec(intvecvec &vec, intvec in){
  vec->vec[vec->veclen++] = in;
}

void pushback_orbvec(orbvec &vec, OrbitMetric_raw in){
  vec->vec[vec->veclen++] = in;
}

void popback_edgevec(edgevec &vec){
  vec->vec[--vec->veclen++] = 0;
}

void popback_intvec(intvec &vec){
  vec->vec[--vec->veclen] = 0;
}

void popback_adjlist(A_Network_raw &vec){
  vec->nodes[--vec->veclen] = 0;
}

void popback_intvecvec(intvecvec &vec){
  vec->vec[--vec->veclen] = 0;
}

void popback_orbvec(orbvec &vec){
  vec->vec[--vec->veclen] = 0;
}

void clear_edgevec(edgevec &vec){
  vec->veclen = 0;
}

void clear_intvec(intvec &vec){
  vec->veclen = 0;
}

void clear_adjlist(A_Network_raw &vec){
  vec->veclen = 0;
}

void clear_intvecvec(intvecvec &vec){
  vec->veclen = 0;
}

void clear_orbvec(orbvec &vec){
  vec->veclen = 0;
}

#endif
