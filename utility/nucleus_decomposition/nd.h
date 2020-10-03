#include <iostream>
#include <sstream>
#include <fstream>
#include <queue>
#include <stack>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <string>
#include <initializer_list>

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <random>
#include <chrono>
#include <sys/stat.h>
#include "bucket.h"

using namespace std;
typedef long long lol;
typedef pair<int, int> vp;
typedef tuple<int, int, int> vt;
typedef vector<vector<int> > nd_Graph;

struct nd_tree_node;

struct subcore {
	bool visible;
	int rank;
	int K;
	int parent;
	int root;
	vector<int> children;
	int size;
	int nEdge;
	double ed;

	subcore (int k) {
		K = k;
		rank = 0;
		parent = -1;
		root = -1;
		visible = true;
		size = 0;
		nEdge = 0;
		ed = -1;
	}
};

struct helpers {
	helpers (vector<vp>* ael) {
		el = ael;
	}
	helpers (vector<vt>* atris) {
		tris = atris;
	}
	helpers () {}

	vector<vp>* el;
	vector<vt>* tris;
};


template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
	std::hash<T> hasher;
	seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
template<typename S, typename T> struct hash<pair<S, T>>
{
	inline size_t operator()(const pair<S, T> & v) const
	{
		size_t seed = 0;
		::hash_combine(seed, v.first);
		::hash_combine(seed, v.second);
		return seed;
	}
};
}




inline bool less_than (int u, int v, nd_Graph& graph) {
	return (graph[u].size() < graph[v].size() || (graph[u].size() == graph[v].size() && u < v));
}

inline void createOrderedIndexEdges (nd_Graph& graph, vector<vp>& el, vector<int>& xel, nd_Graph& orderedGraph) {
	xel.push_back(0);
	orderedGraph.resize(graph.size());
	for (size_t u = 0; u < graph.size(); u++) {
		for (size_t j = 0; j < graph[u].size(); j++) {
			int v = graph[u][j];
			if (less_than (u, v, graph)) {
				orderedGraph[u].push_back(v);
				vp c (u, v);
				el.push_back(c);
			}
		}
		xel.push_back(el.size());
	}
}

inline int getEdgeId (int u, int v, vector<int>& xel, vector<vp>& el, nd_Graph& graph) {

	int a = u, b = v;
	if (less_than (b, a, graph))
		swap (a, b);

	for (int i = xel[a]; i < xel[a+1]; i++)
		if (el[i].second == b)
			return i;

	printf ("getEdgeId returns -1\n");
	exit(1);
}

inline void intersection (vector<int>& a, vector<int>& b, vector<int>& c) {
	size_t i = 0, j = 0;
	while (i < a.size() && j < b.size()) {
		if (a[i] < b[j])
			i++;
		else if (b[j] < a[i])
			j++;
		else {
			c.push_back(a[i]);
			i++;
			j++;
		}
	}
}


void base_kcore(nd_Graph &graph, int nEdge, vector<int> &K, int *maxCore, vector<nd_tree_node> &nd_tree);
void base_ktruss(nd_Graph &graph, int nEdge, vector<int> &K, int *maxtruss, vector<nd_tree_node> &nd_tree);
void base_k34(nd_Graph &graph, int nEdge, vector<int> &K, int *max34, std::vector<nd_tree_node> &nd_tree);

void createSkeleton (int u, initializer_list<int> neighbors, int* nSubcores, vector<int>& K, vector<subcore>& skeleton,	vector<int>& component, vector<int>& unassigned, vector<vp>& relations);
void updateUnassigned (int t, vector<int>& component, int* cid, vector<vp>& relations, vector<int>& unassigned);
void buildHierarchy (int cn, vector<vp>& relations, vector<subcore>& skeleton, int* nSubcores, int nEdge, int nVtx);
void presentNuclei(int variant, vector<subcore> &skeleton, vector<int> &component, nd_Graph &graph, int nEdge,
                   helpers &ax, std::vector<nd_tree_node> &nd_tree);
