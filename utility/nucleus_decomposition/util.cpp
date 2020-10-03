#include "nd_interface.h"

inline int commons (vector<int>& a, vector<int>& b) {
	int i = 0, j = 0;
	int count = 0;
	while (i < a.size() && j < b.size()) {
		if (a[i] < b[j])
			i++;
		else if (b[j] < a[i])
			j++;
		else {
			count++;
			i++;
			j++;
		}
	}
	return count;
}

void rearrange (vector<subcore>& skeleton) { // rearrange children and parents based on visibility

	for (size_t i = 0; i < skeleton.size(); i++)
		skeleton[i].children.clear();
	for (size_t i = 0; i < skeleton.size() - 1; i++) {
		if (skeleton[i].visible) {
			int pr = skeleton[i].parent;
			while (!skeleton[pr].visible)
				pr = skeleton[pr].parent;
			skeleton[i].parent = pr;
			skeleton[pr].children.push_back (i);
		}
	}
}

// find the r-cliques whose component is index, append the vertices in those cliques to the vertices of its all children, sort, compute its density
void reportSubgraph(int variant, int index, vector<int> &component, helpers &ax, vector<subcore> &skeleton,
					nd_Graph &graph, int nEdge, std::vector<nd_tree_node> &nd_tree,
					std::unordered_map<int, int> &skeleton_to_nd_tree, std::vector<bool> visited) {


	if (skeleton[index].parent == -1 && skeleton[index].K == 0) {
		skeleton[index].size = graph.size();
		skeleton[index].nEdge = nEdge;
		skeleton[index].ed = 0;
		return;
	}

	uint32_t r = 0;
	uint32_t s = 0;

	if (variant == 12) {
		r = 1;
		s = 2;
	}
	else if (variant == 23) {
		r = 2;
		s = 3;
	}
	else if (variant == 34) {
		r = 3;
		s = 4;
	}

	std::fill(visited.begin(), visited.end(), false);

	vector<int> vset;
	if (variant == 12 || variant == 13 || variant == 14) {
		for (int i = 0; i < component.size(); i++) {
			if (component[i] == index) {
				if (!visited[i]) {
					vset.push_back(i);
					visited[i] = true;
				}
			}
		}
	}
	else if (variant == 23 || variant == 24) {
		for (int i = 0; i < component.size(); i++) {
			if (component[i] == index) {
				int u = get<0>((*ax.el)[i]);
				int v = get<1>((*ax.el)[i]);

				if (!visited[u]) {
					vset.push_back(u);
					visited[u] = true;
				}

				if (!visited[v]) {
					vset.push_back(v);
					visited[v] = true;
				}
			}
		}
	}
	else if (variant == 34) {
		for (int i = 0; i < component.size(); i++) {
			if (component[i] == index) {
				int u = get<0>((*ax.tris)[i]);
				int v = get<1>((*ax.tris)[i]);
				int w = get<2>((*ax.tris)[i]);
				if (!visited[u]) {
					vset.push_back(u);
					visited[u] = true;
				}

				if (!visited[v]) {
					vset.push_back(v);
					visited[v] = true;
				}

				if (!visited[w]) {
					vset.push_back (w);
					visited[w] = true;
				}
			}
		}
	}

	for (auto child : skeleton[index].children) {
		int nd_node_id = skeleton_to_nd_tree[child];

		for (auto u : nd_tree[nd_node_id].vertices_) {
			if (!visited[u]) {
				vset.push_back(u);
				visited[u] = true;
			}
		}
	}

	std::sort(vset.begin(), vset.end());

	// edge density
	int edge_count = 0;
	for (size_t i = 0; i < vset.size(); i++)
		edge_count += commons (vset, graph[vset[i]]);

	edge_count /= 2;

	double density = 0;
	if (vset.size() > 1)
		density = edge_count /(vset.size() * (vset.size() - 1) / 2.0);

	nd_tree_node node;
	node.parent_ = -1;
	node.k_ = skeleton[index].K;
	node.r_ = r;
	node.s_ = s;
	node.num_edges_ = edge_count;
	node.density_ = density;
	node.vertices_ = vset;
	node.id_ = (int)nd_tree.size();
	skeleton_to_nd_tree[index] = node.id_;
	for (auto child : skeleton[index].children) {
		int node_id = skeleton_to_nd_tree[child];
		node.children_.push_back(node_id);
		nd_tree[node_id].parent_ = node.id_;
	}
	nd_tree.push_back(node);
}

void bfsHierarchy (vector<subcore>& skeleton, stack<int>& scs) {
	rearrange (skeleton);
	queue<int> bfsorder; // we are doing bfs on the hierarchy tree and push the dequeued nodes to the stack
	bfsorder.push(skeleton.size() - 1);
	while (!bfsorder.empty()) {
		int s = bfsorder.front();
		bfsorder.pop();
		scs.push (s);
		for (int r : skeleton[s].children)
			bfsorder.push (r);
	}
}

inline void findRepresentative (int* child, vector<subcore>& skeleton) {
	int u = *child;
	if (skeleton[u].parent != -1) {
		int pr = skeleton[u].parent;
		while (skeleton[u].K == skeleton[pr].K) {
			u = pr;
			if (skeleton[u].parent != -1)
				pr = skeleton[u].parent;
			else
				break;
		}
	}
	*child = u;
}

void presentNuclei(int variant, vector<subcore> &skeleton, vector<int> &component, nd_Graph &graph, int nEdge,
                   helpers &ax, std::vector<nd_tree_node> &nd_tree) {
	// assign unassigned items to top subcore
	for (int i = 0; i < component.size(); i++)
		if (component[i] == -1)
			component[i] = skeleton.size() - 1;

	// match each component with its representative
	unordered_map<int, int> replace;
	for (int i = 0; i < skeleton.size(); i++) {
		int sc = i;
		int original = sc;
		findRepresentative (&sc, skeleton);
		if (original != sc)
			skeleton[original].visible = false;
		replace[original] = sc;
	}

	// each component takes its representative's component number
	for (int i = 0; i < component.size(); i++)
		if (replace.find (component[i]) != replace.end())
			component[i] = replace[component[i]];


    // Rebuild the tree by skipping the invisible nodes and visit the tree by the bottom-up order.
    stack<int> subcoreStack;
	bfsHierarchy (skeleton, subcoreStack);

	std::vector<bool> visited(graph.size());
	std::unordered_map<int, int> skeleton_to_nd_tree;
	while (!subcoreStack.empty()) {
		int i = subcoreStack.top();
		subcoreStack.pop();
		if (skeleton[i].visible) {
			reportSubgraph(variant, i, component, ax, skeleton, graph, nEdge, nd_tree,
					skeleton_to_nd_tree, visited);
		}
	}
}
