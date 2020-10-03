#include "nd.h"

inline void increment (int u, int v, int w, vector<int>& xtris, vector<vp>& el, vector<int>& xel, nd_Graph& orderedTris, nd_Graph& graph, nd_Graph& FC) {
	int i = -1;
	for (i = xel[u]; i < xel[u+1]; i++)
		if (el[i].second == v)
			break;
	bool fl = false;
	for (int j = 0; j < orderedTris[i].size(); j++)
		if (orderedTris[i][j] == w) {
			FC[i][j]++;
			fl = true;
			break;
		}
}

inline void shake (int m, int u, int v, int w, int* a, int* b, nd_Graph& graph) {
	if (m == u) {
		*a = v;
		*b = w;
		if (less_than (w, v, graph))
			swap (*a, *b);
	}
}

inline int maxOrdered (int a, int b, int c, nd_Graph& graph) {
	if (less_than (b, a, graph))
		return less_than (a, c, graph) ? c : a;
	else
		return less_than (b, c, graph) ? c : b;
}

inline void threeWay (vector<int>& x, vector<int>& y, vector<int>& z, vector<int>& commonNeighbors) {
	int i = 0, j = 0, k = 0;
	while (i < x.size() && j < y.size() && k < z.size()) {
		int a = x[i];
		int b = y[j];
		int c = z[k];
		if (a == b && a == c) {
			commonNeighbors.push_back(a);
			i++; j++; k++;
		}
		else {
			int m = max ({a, b, c});
			if (a != m)
				i++;
			if (b != m)
				j++;
			if (c != m)
				k++;
		}
	}
}

inline int getTriangleId (int u, int v, int w, vector<int>& xtris, vector<vp>& el, vector<int>& xel, nd_Graph& orderedTris, nd_Graph& graph) {

	// given the neighbor triangle u, v, w; get el-id by smallest and middle, then get tris id by largest's index in orderedTris[el-id] : xtris[el-id]+index
	int a, b, m = maxOrdered (u, v, w, graph);
	shake (m, u, v, w, &a, &b, graph);
	shake (m, v, u, w, &a, &b, graph);
	shake (m, w, u, v, &a, &b, graph);

	int i = -1;
	for (i = xel[a]; i < xel[a+1]; i++)
		if (el[i].second == b)
			break;

	if (i == -1) {
		printf ("BUG in el of getTriangleId\n");
		exit(1);
	}

	bool fl = false;
	for (int j = 0; j < orderedTris[i].size(); j++)
		if (orderedTris[i][j] == m)
			return xtris[i] + j;

	if (!fl) {
		printf ("BUG in tris of getTriangleId\n");
		exit(1);
	}
}

void create_triangleList (nd_Graph& orderedGraph, vector<vp>& el, nd_Graph& orderedTris, vector<vt>& tris, vector<int>& xtris, nd_Graph& FC) {

	xtris.push_back(0);
	for (size_t i = 0; i < el.size(); i++) {
		int u = get<0>(el[i]);
		int v = get<1>(el[i]);
		vector<int> commonNeighbors;
		intersection (orderedGraph[u], orderedGraph[v], commonNeighbors);
		for (auto w : commonNeighbors) {
			orderedTris[i].push_back (w);
			FC[i].push_back (0);
			vt tr = make_tuple (u, v, w);
			tris.push_back (tr);
		}
		xtris.push_back(tris.size());
	}
}

// per triangle
lol count4cliques (nd_Graph& graph, nd_Graph& orderedGraph, vector<vp>& el, vector<int>& xel, nd_Graph& orderedTris, vector<vt>& tris, vector<int>& xtris, nd_Graph& FC) {

	lol fc = 0;
	for (auto t : tris) {
		int i = 0, j = 0, k = 0;
		int u = get<0>(t);
		int v = get<1>(t);
		int w = get<2>(t);

		while (i < orderedGraph[u].size() && j < orderedGraph[v].size() && k < orderedGraph[w].size()) {
			int a = orderedGraph[u][i];
			int b = orderedGraph[v][j];
			int c = orderedGraph[w][k];

			if (a == b && a == c) {
				int x = a;
				increment (u, v, w, xtris, el, xel, orderedTris, graph, FC);
				increment (u, v, x, xtris, el, xel, orderedTris, graph, FC);
				increment (u, w, x, xtris, el, xel, orderedTris, graph, FC);
				increment (v, w, x, xtris, el, xel, orderedTris, graph, FC);
				i++; j++; k++;
				fc++;
			}
			else {
				int m = max ({a, b, c});
				if (a != m)
					i++;
				if (b != m)
					j++;
				if (c != m)
					k++;
			}
		}
	}
	return fc;
}

void base_k34(nd_Graph &graph, int nEdge, vector<int> &K, int *max34, std::vector<nd_tree_node> &nd_tree) {

	int nVtx = graph.size();

	// Create directed graph from low degree vertices to higher degree vertices AND prepare a CSR-like structure to index the edges
	vector<vp> el;
	vector<int> xel;
	nd_Graph orderedGraph;
	createOrderedIndexEdges (graph, el, xel, orderedGraph);

	// Creating triangle list
	vector<vector<int>> orderedTris (nEdge); // like orderedGraph: each vector<int> is the list of third vertices, w, in the triangles of i-th edge, u - v, in el s.t. u < v < w (deg ordering)
	vector<vt> tris; // like el: list of the triangles aligned to the order in orderedTriangles
	vector<int> xtris; // like xel: indices in tris that starts the triangle list for an edge
	vector<vector<int>> FC (nEdge); // 4-clique counts of each triangle in the orderedTris structure
	create_triangleList (orderedGraph, el, orderedTris, tris, xtris, FC);

	// 4-clique counting for each triangle
	lol fc = count4cliques (graph, orderedGraph, el, xel, orderedTris, tris, xtris, FC);

	// Peeling
	K.resize (tris.size(), -1);
	Naive_Bucket nBucket;
	nBucket.Initialize (nEdge, tris.size()); // maximum 4-clique count of a triangle is nVtx

	int id = 0;
	for (size_t i = 0; i < FC.size(); i++)
		for (size_t j = 0; j < FC[i].size(); j++) {
			if (FC[i][j] > 0)
				nBucket.Insert (xtris[i]+j, FC[i][j]);
			else
				K[xtris[i] + j] = 0;
			id++;
		}

	int fc_t = 0;

	// required for hierarchy
	int cid; // subcore id number
	vector<subcore> skeleton; // equal K valued cores
	vector<int> component; // subcore ids for each int
	vector<vp> relations;
	vector<int> unassigned;
	int nSubcores;


	cid = 0;
	nSubcores = 0;
	component.resize (tris.size(), -1);

	while (true) {
		int t;
		int val ;
		if (nBucket.PopMin(&t, &val)) // if the bucket is empty
			break;

		unassigned.clear();
		subcore sc (val);
		skeleton.push_back (sc);


		fc_t = K[t] = val;

		int u = get<0> (tris[t]);
		int v = get<1> (tris[t]);
		int w = get<2> (tris[t]);
		vector<int> commonNeighbors;
		threeWay (graph[u], graph[v], graph[w], commonNeighbors);
		for (auto x : commonNeighbors) { // decrease the FC of the neighbor triangles with greater FC
			int p = getTriangleId (u, v, x, xtris, el, xel, orderedTris, graph);
			int r = getTriangleId (u, w, x, xtris, el, xel, orderedTris, graph);
			int s = getTriangleId (v, w, x, xtris, el, xel, orderedTris, graph);
			if (K[p] == -1 && K[r] == -1 && K[s] == -1) {
				if (nBucket.CurrentValue(p) > fc_t)
					nBucket.DecVal(p);
				if (nBucket.CurrentValue(r) > fc_t)
					nBucket.DecVal(r);
				if (nBucket.CurrentValue(s) > fc_t)
					nBucket.DecVal(s);
			}
			else
				createSkeleton (t, {p, r, s}, &nSubcores, K, skeleton, component, unassigned, relations);
		}

		updateUnassigned (t, component, &cid, relations, unassigned);
	}

	nBucket.Free();
	*max34 = fc_t; // fc_t is fc of the last popped triangle

	const auto p2 = chrono::steady_clock::now();


	buildHierarchy (*max34, relations, skeleton, &nSubcores, nEdge, nVtx);

	helpers hp (&tris);
	presentNuclei(34, skeleton, component, graph, nEdge, hp, nd_tree);

}
