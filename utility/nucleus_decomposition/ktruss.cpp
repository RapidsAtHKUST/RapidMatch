#include "nd.h"

inline int checkConnectedness (nd_Graph& graph, nd_Graph& orderedGraph, nd_Graph& TC, int u, int v, vector<int>* xel = NULL) {
	int a = u, b = v;
	if (less_than (b, a, graph))
		swap (a, b);
	for (size_t k = 0; k < orderedGraph[a].size(); k++)
		if (orderedGraph[a][k] == b) {
			TC[a][k]++;
			if (xel == NULL)
				return b;
			else
				return (*xel)[a] + k;
		}
	return -1;
}

// per edge
lol countTriangles (nd_Graph& graph, nd_Graph& orderedGraph, nd_Graph& TC) {
    lol tc = 0;

	for (size_t i = 0; i < orderedGraph.size(); i++) {
		for (size_t j = 0; j < orderedGraph[i].size(); j++) {
			for (size_t k = j + 1; k < orderedGraph[i].size(); k++) {
				int a = orderedGraph[i][j];
				int b = orderedGraph[i][k];
				int c = checkConnectedness (graph, orderedGraph, TC, a, b);
				if (c != -1) {
					TC[i][j]++;
					TC[i][k]++;
					tc++;
				}
			}
		}
	}
	return tc;
}

void base_ktruss(nd_Graph &graph, int nEdge, vector<int> &K, int *maxtruss, vector<nd_tree_node> &nd_tree) {

	const auto t1 = chrono::steady_clock::now();
	int nVtx = graph.size();

	// Create directed graph from low degree vertices to higher degree vertices AND prepare a CSR-like structure to index the edges
	vector<vp> el;
	vector<int> xel;
	nd_Graph orderedGraph;
	createOrderedIndexEdges (graph, el, xel, orderedGraph);

	// Triangle counting for each edge
	vector<vector<int> > TC (nVtx);
	for (int i = 0; i < nVtx; i++)
		TC[i].resize (orderedGraph[i].size(), 0);

    countTriangles (graph, orderedGraph, TC);

	// Peeling
	K.resize (nEdge, -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (nVtx, nEdge); // maximum triangle count of an edge is nVtx
	int id = 0;
	for (size_t i = 0; i < orderedGraph.size(); i++)
		for (size_t j = 0; j < orderedGraph[i].size(); j++) {
			if (TC[i][j] > 0)
				nBucket.Insert (id++, TC[i][j]);
			else
				K[id++] = 0;
		}

	int tc_e = 0;

	// required for hierarchy
	int cid; // subcore id number
	vector<subcore> skeleton; // equal K valued cores
	vector<int> component; // subcore ids for each int
	vector<vp> relations;
	vector<int> unassigned;
	int nSubcores;

	cid = 0;
	nSubcores = 0;
	component.resize (nEdge, -1);


	while (true) {
		int e;
		int val;
		if (nBucket.PopMin(&e, &val) == -1) // if the bucket is empty
			break;

		unassigned.clear();
		subcore sc (val);
		skeleton.push_back (sc);


		tc_e = K[e] = val;

		int u = el[e].first;
		int v = el[e].second;
		vector<int> commonNeighbors;
		intersection (graph[u], graph[v], commonNeighbors);
		for (auto w : commonNeighbors) { // decrease the TC of the neighbor edges with greater TC
			int f = getEdgeId (u, w, xel, el, graph);
			int g = getEdgeId (v, w, xel, el, graph);
			if (K[f] == -1 && K[g] == -1) {
				if (nBucket.CurrentValue(f) > tc_e)
					nBucket.DecVal(f);
				if (nBucket.CurrentValue(g) > tc_e)
					nBucket.DecVal(g);
			}
			else
				createSkeleton (e, {f, g}, &nSubcores, K, skeleton, component, unassigned, relations);
		}


		updateUnassigned (e, component, &cid, relations, unassigned);
	}

	nBucket.Free();
	*maxtruss = tc_e;

	buildHierarchy (*maxtruss, relations, skeleton, &nSubcores, nEdge, nVtx);
	helpers hp (&el);
	presentNuclei(23, skeleton, component, graph, nEdge, hp, nd_tree);
}