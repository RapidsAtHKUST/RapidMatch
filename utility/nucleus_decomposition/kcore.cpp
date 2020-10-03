
#include "nd_interface.h"

void base_kcore(nd_Graph &graph, int nEdge, vector<int> &K, int *maxCore, std::vector<nd_tree_node> &nd_tree) {

	const auto p1 = chrono::steady_clock::now();
	int nVtx = graph.size();
	int maxDeg = 0;
	for (auto g : graph)
		if (g.size() > maxDeg)
			maxDeg = g.size();

	// Peeling
	K.resize (graph.size(), -1);
	Naive_Bucket nBucket;
	nBucket.Initialize(maxDeg, nVtx);
	for (int i = 0; i < graph.size(); i++) {
		if (graph[i].size() > 0)
			nBucket.Insert (i, graph[i].size());
		else
			K[i] = 0;
	}

	int deg_u = 0;

	// required for hierarchy
	int cid; // subcore id number
	vector<subcore> skeleton; // equal K valued cores
	vector<int> component; // subcore ids for each int
	vector<vp> relations;
	vector<int> unassigned;
	int nSubcores;


	cid = 0;
	nSubcores = 0;
	component.resize (graph.size(), -1);


	while (true) {
		int u, val;
		if (nBucket.PopMin(&u, &val) == -1) // if the bucket is empty
			break;


		unassigned.clear();
		subcore sc (val);
		skeleton.push_back (sc);


		deg_u = K[u] = val;

		for (auto v : graph[u]) { // decrease the degree of the neighbors with greater degree
			int deg_v = nBucket.CurrentValue(v);
			if (deg_v > deg_u)
				nBucket.DecVal(v);
			else
				createSkeleton (u, {v}, &nSubcores, K, skeleton, component, unassigned, relations);
		}

		updateUnassigned (u, component, &cid, relations, unassigned);
	}

	nBucket.Free();
	*maxCore = deg_u; // deg_u is degree of the last popped int

	buildHierarchy (*maxCore, relations, skeleton, &nSubcores, nEdge, nVtx);
	helpers hp;
	presentNuclei(12, skeleton, component, graph, nEdge, hp, nd_tree);
}
