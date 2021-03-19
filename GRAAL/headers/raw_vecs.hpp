#ifndef RAW_VECS_HPP
#define RAW_VECS_HPP 

#include <stdlib.h>

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
};

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

intvec new_intvec(int cap) {
    intvec new_vec;
    new_vec.vec = (int*) malloc(cap*sizeof(int));
    new_vec.veclen = 0;
    return new_vec;
}

edgevec new_edgevec(int cap) {
    edgevec new_vec;
    new_vec.vec = (Edge*) malloc(cap*sizeof(struct Edge));
    new_vec.veclen = 0;
    return new_vec;
}

Adjlist new_adjlist(int list_cap, int ops_cap) {
    Adjlist new_adjlist;
    new_adjlist.Row = 0;
    new_adjlist.list_len = 0;
    new_adjlist.ops_len = 0;
    new_adjlist.ListW = new_edgevec(list_cap);
    new_adjlist.Ops = new_intvec(ops_cap);
    return new_adjlist;
}

intvecvec new_intvecvec(int cap) {
    intvecvec new_vec;
    new_vec.vec = (intvec*) malloc(cap*sizeof(struct intvec));
    new_vec.veclen = 0;
    return new_vec;
}

A_Network_raw new_network(int cap) {
    A_Network_raw new_net;
    new_net.vec = (Adjlist*) malloc(cap*sizeof(struct Adjlist));
    new_net.nodes_len = 0;
    return new_net;
}

orbvec new_orbvec(int cap) {
    orbvec new_vec;
    new_vec.vec = (OrbitMetric_raw*) malloc(cap*sizeof(struct OrbitMetric_raw));
    new_vec.veclen = 0;
    return new_vec;
}

void delete_intvec(intvec& vec) {
    free(vec.vec);
    vec.vec = NULL;
    vec.veclen = 0;
}

void delete_edgevec(edgevec& vec) {
    free(vec.vec);
    vec.vec = NULL;
    vec.veclen = 0;
}

void delete_adjlist(Adjlist& adj) {
    delete_edgevec(adj.ListW);
    delete_intvec(adj.Ops);
    adj.Row = 0;
    adj.list_len = 0;
    adj.ops_len = 0;
}

void delete_intvecvec(intvecvec& vec) {
    free(vec.vec);
    vec.vec = NULL;
    vec.veclen = 0;
}

void delete_network(A_Network_raw& net) {
    free(net.vec);
    net.vec = NULL;
    net.nodes_len = 0;
}

void delete_orbvec(orbvec& vec) {
    free(vec.vec);
    vec.vec = NULL;
    vec.veclen = 0;
}

void pushback_edgevec(edgevec &vec, Edge in){
  vec.vec[vec.veclen++] = in;
}

void pushback_intvec(intvec &vec, int in){
  vec.vec[vec.veclen++] = in;
}

void pushback_adjlist(A_Network_raw &vec, Adjlist in){
  vec.vec[vec.nodes_len++] = in;
}

void pushback_intvecvec(intvecvec &vec, intvec in){
  vec.vec[vec.veclen++] = in;
}

void pushback_orbvec(orbvec &vec, OrbitMetric_raw in){
  vec.vec[vec.veclen++] = in;
}

void clear_edgevec(edgevec &vec){
  vec.veclen = 0;
}

void clear_intvec(intvec &vec){
  vec.veclen = 0;
}

void clear_adjlist(A_Network_raw &vec){
  vec.nodes_len = 0;
}

void clear_intvecvec(intvecvec &vec){
  vec.veclen = 0;
}

void clear_orbvec(orbvec &vec){
  vec.veclen = 0;
}

void popback_edgevec(edgevec &vec){
  vec.vec[--vec.veclen].first = 0;
  vec.vec[vec.veclen].second = 0;
}

void popback_intvec(intvec &vec){
  vec.vec[--vec.veclen] = 0;
}

void popback_adjlist(A_Network_raw &vec){
  vec.vec[--vec.nodes_len].Row = 0;
  vec.vec[vec.nodes_len].list_len = 0;
  vec.vec[vec.nodes_len].ops_len = 0;
  clear_edgevec(vec.vec[vec.nodes_len].ListW);
  clear_intvec(vec.vec[vec.nodes_len].Ops);
}

void popback_intvecvec(intvecvec &vec){
  clear_intvec(vec.vec[--vec.veclen]);
}

void popback_orbvec(orbvec &vec){
  vec.vec[--vec.veclen].orbitNumber = 0;
  clear_intvec(vec.vec[vec.veclen].orbitDegree);
  clear_intvec(vec.vec[vec.veclen].orbitDistance);
}

#endif
