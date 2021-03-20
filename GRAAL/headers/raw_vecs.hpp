#ifndef RAW_VECS_HPP
#define RAW_VECS_HPP 

#include <umpire/Allocator.hpp>

#include <stdlib.h>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define FIDO_HOST_DEVICE __host__ __device__
#define FIDO_DEVICE __device__
#define FIDO_HOST __host__
#define FIDO_CONSTANT __constant__
#else
#define FIDO_HOST_DEVICE
#define FIDO_DEVICE
#define FIDO_HOST
#define FIDO_CONSTANT
#endif

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

FIDO_HOST_DEVICE intvec new_intvec(int cap) {
	intvec new_vec;
	new_vec.vec = (int*) malloc(cap*sizeof(int));
	new_vec.veclen = 0;
	return new_vec;
}

intvec new_intvec_umpire(int cap, umpire::Allocator& alloc) {
	intvec new_vec;
	new_vec.vec = static_cast<int*>(alloc.allocate(cap*sizeof(int)));
	new_vec.veclen = 0;
	return new_vec;
}

FIDO_HOST_DEVICE edgevec new_edgevec(int cap) {
	edgevec new_vec;
	new_vec.vec = (Edge*) malloc(cap*sizeof(struct Edge));
	new_vec.veclen = 0;
	return new_vec;
}

edgevec new_edgevec_umpire(int cap, umpire::Allocator& alloc) {
	edgevec new_vec;
	new_vec.vec = static_cast<Edge*>(alloc.allocate(cap*sizeof(struct Edge)));
	new_vec.veclen = 0;
	return new_vec;
}

FIDO_HOST_DEVICE Adjlist new_adjlist(int list_cap, int ops_cap) {
	Adjlist new_adjlist;
	new_adjlist.Row = 0;
	new_adjlist.list_len = 0;
	new_adjlist.ops_len = 0;
	new_adjlist.ListW = new_edgevec(list_cap);
	new_adjlist.Ops = new_intvec(ops_cap);
	return new_adjlist;
}

Adjlist new_adjlist_umpire(int list_cap, int ops_cap, umpire::Allocator& alloc) {
	Adjlist new_adjlist;
	new_adjlist.Row = 0;
	new_adjlist.list_len = 0;
	new_adjlist.ops_len = 0;
	new_adjlist.ListW = new_edgevec_umpire(list_cap, alloc);
	new_adjlist.Ops = new_intvec_umpire(ops_cap, alloc);
	return new_adjlist;
}

FIDO_HOST_DEVICE intvecvec new_intvecvec(int cap) {
	intvecvec new_vec;
	new_vec.vec = (intvec*) malloc(cap*sizeof(struct intvec));
	new_vec.veclen = 0;
	return new_vec;
}

intvecvec new_intvecvec_umpire(int cap, umpire::Allocator& alloc) {
	intvecvec new_vec;
	new_vec.vec = static_cast<intvec*>(alloc.allocate(cap*sizeof(struct intvec)));
	new_vec.veclen = 0;
	return new_vec;
}

FIDO_HOST_DEVICE A_Network_raw new_network(int cap) {
	A_Network_raw new_net;
	new_net.vec = (Adjlist*) malloc(cap*sizeof(struct Adjlist));
	new_net.nodes_len = 0;
	return new_net;
}

A_Network_raw new_network_umpire(int cap, umpire::Allocator& alloc) {
	A_Network_raw new_net;
	new_net.vec = static_cast<Adjlist*>(alloc.allocate(cap*sizeof(struct Adjlist)));
	new_net.nodes_len = 0;
	return new_net;
}

FIDO_HOST_DEVICE orbvec new_orbvec(int cap) {
	orbvec new_vec;
	new_vec.vec = (OrbitMetric_raw*) malloc(cap*sizeof(struct OrbitMetric_raw));
	new_vec.veclen = 0;
	return new_vec;
}

orbvec new_orbvec_umpire(int cap, umpire::Allocator& alloc) {
	orbvec new_vec;
	new_vec.vec = static_cast<OrbitMetric_raw*>(alloc.allocate(cap*sizeof(struct OrbitMetric_raw)));
	new_vec.veclen = 0;
	return new_vec;
}

FIDO_HOST_DEVICE void delete_intvec(intvec& vec) {
	free(vec.vec);
	vec.vec = NULL;
	vec.veclen = 0;
}

void delete_intvec_umpire(intvec& vec, umpire::Allocator& alloc) {
	alloc.deallocate(vec.vec);
	vec.vec = NULL;
	vec.veclen = 0;
}

FIDO_HOST_DEVICE void delete_edgevec(edgevec& vec) {
	free(vec.vec);
	vec.vec = NULL;
	vec.veclen = 0;
}

void delete_edgevec_umpire(edgevec& vec, umpire::Allocator& alloc) {
	alloc.deallocate(vec.vec);
	vec.vec = NULL;
	vec.veclen = 0;
}

FIDO_HOST_DEVICE void delete_adjlist(Adjlist& adj) {
	delete_edgevec(adj.ListW);
	delete_intvec(adj.Ops);
	adj.Row = 0;
	adj.list_len = 0;
	adj.ops_len = 0;
}

void delete_adjlist_umpire(Adjlist& adj, umpire::Allocator& alloc) {
	delete_edgevec_umpire(adj.ListW, alloc);
	delete_intvec_umpire(adj.Ops, alloc);
	adj.Row = 0;
	adj.list_len = 0;
	adj.ops_len = 0;
}

FIDO_HOST_DEVICE void delete_intvecvec(intvecvec& vec) {
	free(vec.vec);
	vec.vec = NULL;
	vec.veclen = 0;
}

void delete_intvecvec_umpire(intvecvec& vec, umpire::Allocator& alloc) {
	alloc.deallocate(vec.vec);
	vec.vec = NULL;
	vec.veclen = 0;
}

FIDO_HOST_DEVICE void delete_network(A_Network_raw& net) {
	free(net.vec);
	net.vec = NULL;
	net.nodes_len = 0;
}

void delete_network_umpire(A_Network_raw& net, umpire::Allocator& alloc) {
	alloc.deallocate(net.vec);
	net.vec = NULL;
	net.nodes_len = 0;
}

FIDO_HOST_DEVICE void delete_orbitmetric(OrbitMetric_raw& orb) {
	delete_intvec(orb.orbitDegree);
	delete_intvec(orb.orbitDistance);
	orb.orbitNumber = 0;
}

void delete_orbitmetric_umpire(OrbitMetric_raw& orb, umpire::Allocator& alloc) {
	delete_intvec_umpire(orb.orbitDegree, alloc);
	delete_intvec_umpire(orb.orbitDistance, alloc);
	orb.orbitNumber = 0;
}

FIDO_HOST_DEVICE void delete_orbvec(orbvec& vec) {
	free(vec.vec);
	vec.vec = NULL;
	vec.veclen = 0;
}

void delete_orbvec_umpire(orbvec& vec, umpire::Allocator& alloc) {
	alloc.deallocate(vec.vec);
	vec.vec = NULL;
	vec.veclen = 0;
}

FIDO_HOST_DEVICE void pushback_edgevec(edgevec &vec, Edge in){
	vec.vec[vec.veclen++] = in;
}

FIDO_HOST_DEVICE void pushback_intvec(intvec &vec, int in){
	vec.vec[vec.veclen++] = in;
}

FIDO_HOST_DEVICE void pushback_adjlist(A_Network_raw &vec, Adjlist in){
	vec.vec[vec.nodes_len++] = in;
}

FIDO_HOST_DEVICE void pushback_intvecvec(intvecvec &vec, intvec in){
	vec.vec[vec.veclen++] = in;
}

FIDO_HOST_DEVICE void pushback_orbvec(orbvec &vec, OrbitMetric_raw in){
	vec.vec[vec.veclen++] = in;
}

FIDO_HOST_DEVICE void clear_edgevec(edgevec &vec){
	vec.veclen = 0;
}

FIDO_HOST_DEVICE void clear_intvec(intvec &vec){
	vec.veclen = 0;
}

FIDO_HOST_DEVICE void clear_adjlist(A_Network_raw &vec){
	for (int i = 0; i < vec.nodes_len; i++)
	{
		delete_adjlist(vec.vec[i]);
	}
	vec.nodes_len = 0;
}

FIDO_HOST_DEVICE void clear_intvecvec(intvecvec &vec){
	for (int i = 0; i < vec.veclen; i++)
	{
		delete_intvec(vec.vec[i]);
	}
	vec.veclen = 0;
}

FIDO_HOST_DEVICE void clear_orbvec(orbvec &vec){
	for (int i = 0; i < vec.veclen; i++)
	{
		delete_orbitmetric(vec.vec[i]);
	}
	vec.veclen = 0;
}

FIDO_HOST_DEVICE void popback_edgevec(edgevec &vec){
	vec.vec[--vec.veclen].first = 0;
	vec.vec[vec.veclen].second = 0;
}

FIDO_HOST_DEVICE void popback_intvec(intvec &vec){
	vec.vec[--vec.veclen] = 0;
}

FIDO_HOST_DEVICE void popback_adjlist(A_Network_raw &vec){
	vec.vec[--vec.nodes_len].Row = 0;
	vec.vec[vec.nodes_len].list_len = 0;
	vec.vec[vec.nodes_len].ops_len = 0;
	clear_edgevec(vec.vec[vec.nodes_len].ListW);
	clear_intvec(vec.vec[vec.nodes_len].Ops);
}

FIDO_HOST_DEVICE void popback_intvecvec(intvecvec &vec){
	clear_intvec(vec.vec[--vec.veclen]);
}

FIDO_HOST_DEVICE void popback_orbvec(orbvec &vec){
	vec.vec[--vec.veclen].orbitNumber = 0;
	clear_intvec(vec.vec[vec.veclen].orbitDegree);
	clear_intvec(vec.vec[vec.veclen].orbitDistance);
}

FIDO_HOST_DEVICE void copy_intvec(intvec& src, intvec& dst) {
	for (int i = 0; i < src.veclen; i++)
	{
		dst.vec[i] = src.vec[i];
	}
	dst.veclen = src.veclen;
}

FIDO_HOST_DEVICE void copy_intvecvec(intvecvec& src, intvecvec& dst) {
	for (int i = 0; i < src.veclen; i++)
	{
		copy_intvec(src.vec[i], dst.vec[i]);
	}
	dst.veclen = src.veclen;
}

FIDO_HOST_DEVICE void copy_edgevec(edgevec& src, edgevec& dst) {
	for (int i = 0; i < src.veclen; i++)
	{
		dst.vec[i] = src.vec[i];
	}
	dst.veclen = src.veclen;
}

void stdnet_to_rawnet(A_Network &stdnet, A_Network_raw &rawnet, umpire::Allocator &alloc){
	rawnet = new_network_umpire(stdnet.size(), alloc);
	for(ADJ_Bundle &node : stdnet){
		pushback_adjlist(rawnet, new_adjlist_umpire(node.ListW.size(),0, alloc));
		rawnet.vec[rawnet.nodes_len-1].Row = node.Row;
		memcpy(rawnet.vec[rawnet.nodes_len-1].ListW, node.ListW.data(), node.ListW.size()*sizeof(Edge));
	}
}

void stdorb_to_raworb(vector<OrbitMetric> &stdorb, orbvec &raworb, umpire::Allocator &alloc){
	raworb = new_orbvec_umpire(stdorb.size(), alloc);
	for(OrbitMetric &node : stdorb){
		OrbitMetric_raw tmp;
		tmp.orbitDegree = new_intvec_umpire(node.orbitDegree.size(), alloc);
		tmp.orbitDistance = new_intvec_umpire(node.orbitDistance.size(), alloc);
		tmp.orbitNumber = node.orbitNumber;
		memcpy(tmp.orbitDegree.vec, node.OrbitDegree.data(), node.OrbitDegree.size()*sizeof(int));
		memcpy(tmp.orbitDistance.vec, node.OrbitDistance.data(), node.OrbitDistance.size()*sizeof(int));
		pushback_orbvec(raworb, tmp);
	}
}

#endif
