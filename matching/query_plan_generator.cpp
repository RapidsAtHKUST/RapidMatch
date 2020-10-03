#include "query_plan_generator.h"
#include "utility/computesetintersection.h"
#include <cmath>
#include <bits/stdc++.h>

double query_plan_generator::ordering_time_;
double query_plan_generator::traversal_time_;
double query_plan_generator::nd_time_;

void query_plan_generator::generate_query_plan_for_test(Graph *query_graph, std::vector<uint32_t> &order) {
    uint32_t query_vertices_num = query_graph->getVerticesCount();
    uint32_t core_vertices_num = query_graph->get2CoreSize();
    std::vector<bool> visited(query_vertices_num, false);
    order.resize(query_vertices_num);

    // Select the vertex with the maximum degree as the start vertex.
    order[0] = 0;
    for (uint32_t u = 1; u < query_vertices_num; ++u) {
        if (query_graph->getVertexDegree(u) > query_graph->getVertexDegree(order[0]) && (core_vertices_num == 0 || query_graph->getCoreValue(u)) ) {
            order[0] = u;
        }
    }
    visited[order[0]] = true;
    // Order vertices.
    std::vector<uint32_t > tie_vertices;
    std::vector<uint32_t > temp;

    for (uint32_t  i = 1; i < query_vertices_num; ++i) {
        // Select the vertices with the maximum number of backward neighbors.
        uint32_t  max_bn = 0;
        for (uint32_t  u = 0; u < query_vertices_num; ++u) {
            if (!visited[u] && (query_graph->getCoreValue(u) > 1 || i >= core_vertices_num)) {
                // Compute the number of backward neighbors of u.
                uint32_t  cur_bn = 0;
                for (uint32_t  j = 0; j < i; ++j) {
                    uint32_t  uu = order[j];
                    if (query_graph->checkEdgeExistence(u, uu)) {
                        cur_bn += 1;
                    }
                }

                // Update the vertices under consideration.
                if (cur_bn > max_bn) {
                    max_bn = cur_bn;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                } else if (cur_bn == max_bn) {
                    tie_vertices.push_back(u);
                }
            }
        }

        if (tie_vertices.size() != 1) {
            temp.swap(tie_vertices);
            tie_vertices.clear();

            uint32_t  count = 0;
            std::vector<uint32_t > u_fn;
            for (auto u : temp) {
                // Compute the number of vertices in the matching order that has at least one vertex not in the matching order && connected with u.

                // Get the neighbors of u that are not in the matching order.
                uint32_t  un_count;
                const uint32_t * un = query_graph->getVertexNeighbors(u, un_count);
                for (uint32_t  j = 0; j < un_count; ++j) {
                    if (!visited[un[j]]) {
                        u_fn.push_back(un[j]);
                    }
                }

                // Compute the valid number of vertices.
                uint32_t  cur_count = 0;
                for (uint32_t  j = 0; j < i; ++j) {
                    uint32_t  uu = order[j];
                    uint32_t  uun_count;
                    const uint32_t * uun = query_graph->getVertexNeighbors(uu, uun_count);
                    uint32_t  common_neighbor_count = 0;
                    ComputeSetIntersection::ComputeCandidates(uun, uun_count, u_fn.data(), (uint32_t )u_fn.size(), common_neighbor_count);
                    if (common_neighbor_count > 0) {
                        cur_count += 1;
                    }
                }

                u_fn.clear();

                // Update the vertices under consideration.
                if (cur_count > count) {
                    count = cur_count;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                }
                else if (cur_count == count){
                    tie_vertices.push_back(u);
                }
            }
        }

        if (tie_vertices.size() != 1) {
            temp.swap(tie_vertices);
            tie_vertices.clear();

            uint32_t  count = 0;
            std::vector<uint32_t > u_fn;
            for (auto u : temp) {
                // Compute the number of vertices not in the matching order && not the neighbor of vertices in the matching order, but is connected with u.

                // Get the neighbors of u that are not in the matching order.
                uint32_t  un_count;
                const uint32_t * un = query_graph->getVertexNeighbors(u, un_count);
                for (uint32_t  j = 0; j < un_count; ++j) {
                    if (!visited[un[j]]) {
                        u_fn.push_back(un[j]);
                    }
                }

                // Compute the valid number of vertices.
                uint32_t  cur_count = 0;
                for (auto uu : u_fn) {
                    bool valid = true;

                    for (uint32_t  j = 0; j < i; ++j) {
                        if (query_graph->checkEdgeExistence(uu, order[j])) {
                            valid = false;
                            break;
                        }
                    }

                    if (valid) {
                        cur_count += 1;
                    }
                }

                u_fn.clear();

                // Update the vertices under consideration.
                if (cur_count > count) {
                    count = cur_count;
                    tie_vertices.clear();
                    tie_vertices.push_back(u);
                }
                else if (cur_count == count){
                    tie_vertices.push_back(u);
                }
            }
        }

        order[i] = tie_vertices[0];

        visited[order[i]] = true;
        tie_vertices.clear();
        temp.clear();
    }
}

void query_plan_generator::generate_query_plan_with_nd(Graph *query_graph, catalog *storage,
                                                       std::vector<std::vector<uint32_t>>& vertex_orders) {
    std::vector<nd_tree_node> density_tree;
    std::vector<nd_tree_node> k12_tree;
    std::vector<nd_tree_node> k23_tree;
    std::vector<nd_tree_node> k34_tree;
    auto start = std::chrono::high_resolution_clock::now();

    nd_interface::nd(query_graph, 1, 2, k12_tree);
    nd_interface::nd(query_graph, 2, 3, k23_tree);
    nd_interface::nd(query_graph, 3, 4, k34_tree);

    auto end = std::chrono::high_resolution_clock::now();
    nd_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    construct_density_tree(density_tree, k12_tree, k23_tree, k34_tree);

    std::vector<std::vector<uint32_t>> node_orders;
    traversal_density_tree(query_graph, storage, density_tree, vertex_orders, node_orders);

    end = std::chrono::high_resolution_clock::now();
    traversal_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    ordering_time_ = nd_time_ + traversal_time_;

    nd_interface::print_nd_tree(1, 2, k12_tree);
    nd_interface::print_nd_tree(2, 3, k23_tree);
    nd_interface::print_nd_tree(3, 4, k34_tree);
    print_density_tree(density_tree);
    print_vertex_orders(query_graph, vertex_orders, node_orders);
}

void query_plan_generator::construct_density_tree(std::vector<nd_tree_node> &density_tree,
                                                  std::vector<nd_tree_node> &k12_tree,
                                                  std::vector<nd_tree_node> &k23_tree,
                                                  std::vector<nd_tree_node> &k34_tree) {
    /** Eliminate the duplicate nodes.
     *     1. For each nucleus in k12_tree, if it is a subset of the nucleus in k23, then remove the subtree rooted at it from k12_tree.
     *     2. For each nucleus in k23_tree, if it is a subset of the nucleus in k34, then remove the subtree rooted at it from k23_tree.
     */
    eliminate_node(density_tree, k12_tree, k23_tree);
    eliminate_node(density_tree, k23_tree, k34_tree);

    /**
     * Add the node in k34_tree to the density tree.
     */
    uint32_t offset = density_tree.size();
    for (auto& node : k34_tree) {
        density_tree.push_back(node);
        auto& new_node = density_tree.back();

        if (new_node.parent_ != -1) {
            new_node.parent_ += offset;
        }

        new_node.id_ += offset;

        for (uint32_t i = 0; i < new_node.children_.size(); ++i) {
            new_node.children_[i] += offset;
        }
    }

    /**
     * Construct the tree.
     *     1. For the root node of each k23_tree, find the smallest (by cardinality) 12 nucleus contain it and set the nucleus as its parent.
     *     2. For the root node of each k34_tree, find the smallest (by cardinality) 12 or 23 nucleus contain it and set the nucleus as its parent.
     */
     merge_tree(density_tree);
}

void query_plan_generator::eliminate_node(std::vector<nd_tree_node> &density_tree, std::vector<nd_tree_node> &src_tree,
                                          std::vector<nd_tree_node> &target_tree) {
    std::vector<int> roots;
    for (auto& node : target_tree) {
        if (node.parent_ == -1) {
            roots.push_back(node.id_);
        }
    }

    std::queue<int> q;
    unordered_map<int, int> mappings;
    mappings[-1] = -1;

    for (auto& node : src_tree) {
        if (node.parent_ == -1) {
            q.push(node.id_);
            while (!q.empty()) {
                int cur_node_id = q.front();
                q.pop();

                auto &cur_node = src_tree[cur_node_id];

                bool remove = false;
                for (auto &target_node_id : roots) {
                    auto &target_node = target_tree[target_node_id];
                    if (cur_node.vertices_.size() <= target_node.vertices_.size()) {
                        uint32_t cn_cnt = 0;
                        ComputeSetIntersection::ComputeCandidates((uint32_t *) cur_node.vertices_.data(),
                                                                  (uint32_t) cur_node.vertices_.size(),
                                                                  (uint32_t *) target_node.vertices_.data(),
                                                                  (uint32_t) target_node.vertices_.size(), cn_cnt);

                        if (cn_cnt == cur_node.vertices_.size()) {
                            remove = true;
                        }
                    }
                }

                if (!remove) {
                    auto new_id = static_cast<int>(density_tree.size());
                    mappings[cur_node.id_] = new_id;
                    density_tree.push_back(cur_node);
                    // Update its id, parent id and the children of its parent.
                    density_tree.back().id_ = new_id;
                    density_tree.back().parent_ = mappings[cur_node.parent_];
                    density_tree.back().children_.clear();

                    if (density_tree.back().parent_ != -1)
                        density_tree[density_tree.back().parent_].children_.push_back(new_id);

                    for (auto child_id : cur_node.children_) {
                        q.push(child_id);
                    }
                }
            }
        }
    }
}

void query_plan_generator::merge_tree(std::vector<nd_tree_node> &density_tree) {
    for (int i = density_tree.size() - 1; i >= 0; --i) {
        auto& cur_node = density_tree[i];
        uint32_t smallest_cardinality = std::numeric_limits<uint32_t>::max();
        if (cur_node.parent_ == -1) {
            for (int j = i - 1; j >= 0; --j) {
                auto& target_node = density_tree[j];

                if (target_node.r_ < cur_node.r_ && target_node.vertices_.size() < smallest_cardinality) {
                    uint32_t cn_cnt = 0;
                    ComputeSetIntersection::ComputeCandidates((uint32_t *) cur_node.vertices_.data(),
                                                              (uint32_t) cur_node.vertices_.size(),
                                                              (uint32_t *) target_node.vertices_.data(),
                                                              (uint32_t) target_node.vertices_.size(), cn_cnt);

                    if (cn_cnt == cur_node.vertices_.size()) {
                        cur_node.parent_ = target_node.id_;
                        target_node.children_.push_back(cur_node.id_);
                        smallest_cardinality = target_node.vertices_.size();
                    }
                }

            }
        }
    }
}

void query_plan_generator::print_density_tree(std::vector<nd_tree_node> &density_tree) {
    printf("Density Tree:\n");
    std::queue<int> q;
    q.push(0);

    while (!q.empty()) {
        int node_id = q.front();
        q.pop();
        printf ("\tid %d, parent %d, k %d, r %d, s %d, |V| %zu, |E| %d, density %.4f, num of children %zu: ",
                density_tree[node_id].id_, density_tree[node_id].parent_, density_tree[node_id].k_, density_tree[node_id].r_,
                density_tree[node_id].s_, density_tree[node_id].vertices_.size(), density_tree[node_id].num_edges_,
                density_tree[node_id].density_, density_tree[node_id].children_.size());
        for (auto u : density_tree[node_id].vertices_) {
                printf ("%d ", u);
        }
        printf("\n");

        for (auto child : density_tree[node_id].children_) {
            q.push(child);
        }
    }
}

void query_plan_generator::traversal_density_tree(Graph *query_graph, catalog *storage,
                                                  std::vector<nd_tree_node> &density_tree,
                                                  std::vector<std::vector<uint32_t>> &vertex_orders,
                                                  std::vector<std::vector<uint32_t>> &node_orders) {
    std::vector<uint32_t> vertex_order;
    std::vector<uint32_t> node_order;
    std::vector<bool> visited_vertex(query_graph->getVerticesCount(), false);
    std::vector<bool> visited_node(query_graph->getVerticesCount(), false);
    std::unordered_set<uint32_t> extendable_vertex;

    for (auto& node : density_tree) {
        if (node.children_.empty()) {
            traversal_node(query_graph, storage, density_tree, vertex_order, node_order, visited_vertex, visited_node,
                           node, extendable_vertex);

            vertex_orders.emplace_back(vertex_order);
            node_orders.emplace_back(node_order);

            vertex_order.clear();
            node_order.clear();
            extendable_vertex.clear();

            std::fill(visited_vertex.begin(), visited_vertex.end(), false);
            std::fill(visited_node.begin(), visited_node.end(), false);
        }
    }
}

void query_plan_generator::traversal_node(Graph *query_graph, catalog *storage, std::vector<nd_tree_node> &density_tree,
                                         std::vector<uint32_t> &vertex_order, std::vector<uint32_t> &node_order,
                                         vector<bool> &visited_vertex, std::vector<bool> &visited_node, nd_tree_node &cur_node,
                                         std::unordered_set<uint32_t> &extendable_vertex) {
    node_order.push_back(cur_node.id_);
    visited_node[cur_node.id_] = true;

    // Check if all vertices has been visited.
    bool all_vertices_visited = true;
    for (auto u : cur_node.vertices_) {
        if (!visited_vertex[u]) {
            all_vertices_visited = false;
            break;
        }
    }

    if (all_vertices_visited) {
        if (cur_node.parent_ != -1 && !visited_node[cur_node.parent_])
            traversal_node(query_graph, storage, density_tree, vertex_order, node_order, visited_vertex, visited_node,
                           density_tree[cur_node.parent_], extendable_vertex);
        return;
    }

    std::vector<uint32_t> prev(query_graph->getVerticesCount());
    std::vector<double> dist(query_graph->getVerticesCount());

    // True means that there exists a child that is not visited.
    bool flag = true;
    while (flag) {
        double max_connectivity = 1.5;
        uint32_t selected_child = 0;
        flag = false;

        for (auto child_id : cur_node.children_) {
            if (!visited_node[child_id]) {
                double connectivity = connectivity_common_neighbors(visited_vertex, density_tree[child_id]);
                if (connectivity > max_connectivity
                    || (connectivity == max_connectivity &&
                        density_tree[child_id].density_ > density_tree[selected_child].density_)
                    || (connectivity == max_connectivity &&
                        density_tree[child_id].density_ > density_tree[selected_child].density_ &&
                        density_tree[child_id].vertices_.size() > density_tree[selected_child].vertices_.size())) {
                    max_connectivity = connectivity;
                    selected_child = child_id;
                }

                flag = true;
            }
        }

        if (flag) {
            if (max_connectivity < 2) {
                connectivity_shortest_path(query_graph, storage, vertex_order, visited_vertex, cur_node, prev, dist);
                max_connectivity = 0;
                uint32_t selected_vertex = 0;
                for (auto child_id : cur_node.children_) {
                    if (!visited_node[child_id]) {
                        double local_min_dist = std::numeric_limits<double>::max();
                        uint32_t local_selected_vertex = 0;
                        for (auto u : density_tree[child_id].vertices_) {
                            if (local_min_dist > dist[u]) {
                                local_min_dist = dist[u];
                                local_selected_vertex = u;
                            }
                        }

                        double connectivity = 1.0 / (local_min_dist + 1.0);
                        if (connectivity > max_connectivity) {
                            max_connectivity = connectivity;
                            selected_child = child_id;
                            selected_vertex = local_selected_vertex;
                        }
                    }
                }

                // Construct the path.
                std::vector<uint32_t> vertex_in_shortest_path;
                while (prev[selected_vertex] != selected_vertex) {
                    vertex_in_shortest_path.push_back(selected_vertex);
                    selected_vertex = prev[selected_vertex];
                }

                for (auto rb = vertex_in_shortest_path.rbegin(); rb != vertex_in_shortest_path.rend(); rb++) {
                    if (!visited_vertex[*rb]) {
                        vertex_order.push_back(*rb);
                        visited_vertex[*rb] = true;
                        update_extendable_vertex(query_graph, *rb, extendable_vertex, visited_vertex);

                        {
                            uint32_t bn_cnt_threshold = 2;
                            if (density_tree[selected_child].r_ > bn_cnt_threshold) {
                                bn_cnt_threshold = density_tree[selected_child].r_;
                            }

                            greedy_expand(query_graph, storage, vertex_order, visited_vertex, extendable_vertex,
                                          bn_cnt_threshold);
                        }
                    }
                }
            }

            traversal_node(query_graph, storage, density_tree, vertex_order, node_order, visited_vertex, visited_node,
                           density_tree[selected_child], extendable_vertex);

        }
    }


    if (vertex_order.empty()) {
        uint32_t min_relation_size = std::numeric_limits<uint32_t>::max();
        uint32_t max_degree_sum = 0;
        uint32_t first_vertex = 0;
        uint32_t second_vertex = 0;
        for (auto u : cur_node.vertices_) {
            uint32_t u_deg = query_graph->getVertexDegree(u);
            for (auto v : cur_node.vertices_) {
                uint32_t v_deg = query_graph->getVertexDegree(v);
                uint32_t local_degree_sum = u_deg + v_deg;
                if (u < v && query_graph->checkEdgeExistence(u, v)) {
                    if (local_degree_sum > max_degree_sum ||
                            (local_degree_sum == max_degree_sum &&
                                storage->get_edge_relation_cardinality(u, v) < min_relation_size)) {
                        max_degree_sum = local_degree_sum;
                        min_relation_size = storage->get_edge_relation_cardinality(u, v);
                        first_vertex = u;
                        second_vertex = v;

                        if (u_deg < v_deg) {
                            swap(first_vertex, second_vertex);
                        }
                    }
                }
            }
        }

        vertex_order.push_back(first_vertex);
        visited_vertex[first_vertex] = true;
        update_extendable_vertex(query_graph, first_vertex, extendable_vertex, visited_vertex);

        vertex_order.push_back(second_vertex);
        visited_vertex[second_vertex] = true;
        update_extendable_vertex(query_graph, second_vertex, extendable_vertex, visited_vertex);
    }

    while (!all_vertices_visited) {
        all_vertices_visited = true;

        uint32_t max_bn = 0;
        uint32_t max_degree = 0;
        uint32_t min_candidate_vertex_num = std::numeric_limits<uint32_t>::max();
        uint32_t selected_vertex = 0;
        for (auto u : cur_node.vertices_) {
            if (!visited_vertex[u]) {
                uint32_t nbrs_cnt;
                const uint32_t *nbrs = query_graph->getVertexNeighbors(u, nbrs_cnt);
                uint32_t local_bn_cnt = 0;
                uint32_t u_deg = query_graph->getVertexDegree(u);
                uint32_t local_candidate_vertex_num = storage->num_candidates_[u];

                for (uint32_t i = 0; i < nbrs_cnt; ++i) {
                    uint32_t v = nbrs[i];

                    if (visited_vertex[v]) {
                        local_bn_cnt += 1;
                    }
                }

                if (local_bn_cnt > max_bn ||
                        (local_bn_cnt == max_bn && local_candidate_vertex_num < min_candidate_vertex_num) ||
                    (local_bn_cnt == max_bn &&local_candidate_vertex_num == min_candidate_vertex_num && u_deg > max_degree)) {
                    selected_vertex = u;
                    min_candidate_vertex_num = local_candidate_vertex_num;
                    max_bn = local_bn_cnt;
                    max_degree = u_deg;
                }

                all_vertices_visited = false;
            }
        }

        if (!all_vertices_visited) {
            vertex_order.push_back(selected_vertex);
            visited_vertex[selected_vertex] = true;
            update_extendable_vertex(query_graph, selected_vertex, extendable_vertex, visited_vertex);
        }
    }
    {
        uint32_t bn_cnt_threshold = 2;
        if (cur_node.parent_ != -1 && density_tree[cur_node.parent_].r_ > bn_cnt_threshold) {
            bn_cnt_threshold = density_tree[cur_node.parent_].r_;
        }

        greedy_expand(query_graph, storage, vertex_order, visited_vertex, extendable_vertex, bn_cnt_threshold);
    }
    if (cur_node.parent_ != -1 && !visited_node[cur_node.parent_])
        traversal_node(query_graph, storage, density_tree, vertex_order, node_order, visited_vertex, visited_node,
                       density_tree[cur_node.parent_], extendable_vertex);
}

double query_plan_generator::connectivity_common_neighbors(std::vector<bool> &visited_vertex, nd_tree_node &child_node) {
    double connectivity = 1;
    for (auto u : child_node.vertices_) {
        if (visited_vertex[u]) {
            connectivity += 1;
        }
    }
    return connectivity;
}

void query_plan_generator::connectivity_shortest_path(Graph *query_graph, catalog *storage,
                                                      std::vector<uint32_t> &vertex_order,
                                                      std::vector<bool> &visited_vertex, nd_tree_node &cur_node,
                                                      std::vector<uint32_t> &prev, std::vector<double> &dist) {
    std::fill(prev.begin(), prev.end(), 0);
    std::fill(dist.begin(), dist.end(), std::numeric_limits<double>::max());
    std::vector<bool> in_cur_node(query_graph->getVerticesCount(), false);

    priority_queue<std::pair<double, uint32_t>, std::vector<std::pair<double, uint32_t>>, std::greater<std::pair<double, uint32_t>>> min_pq;

    for (auto u : vertex_order) {
        if (visited_vertex[u]) {
            prev[u] = u;
        }
        min_pq.push(std::make_pair(1, u));
    }

    for (auto u : cur_node.vertices_) {
        in_cur_node[u] = true;
    }

    while (!min_pq.empty()) {
        auto top_element = min_pq.top();
        min_pq.pop();

        uint32_t u = top_element.second;
        double dist_in_pq = top_element.first;

        // If the distance is updated, then skip the element.
        if (dist_in_pq > dist[u])
            continue;

        uint32_t nbrs_cnt;
        const uint32_t* nbrs = query_graph->getVertexNeighbors(u, nbrs_cnt);

        for (uint32_t i = 0; i < nbrs_cnt; ++i) {
            uint32_t v = nbrs[i];
            if (in_cur_node[v] && !visited_vertex[v]) {
                double distance = storage->num_candidates_[v] > 3 ? storage->num_candidates_[v] : 3;
                double updated_dist = dist_in_pq * log(distance);

                if (updated_dist < dist[v]) {
                    dist[v] = updated_dist;
                    prev[v] = u;
                    min_pq.push(std::make_pair(updated_dist, v));
                }
            }
        }
    }
}

void query_plan_generator::print_vertex_orders(Graph *query_graph, std::vector<std::vector<uint32_t>> &vertex_orders,
                                               std::vector<std::vector<uint32_t>> &node_orders) {

    std::vector<std::vector<uint32_t>> selected_vertex_order;
    uint32_t max_utility_value = 0;
    uint32_t selected_vertex_order_id = 0;
    for (uint32_t i = 0; i < vertex_orders.size(); ++i) {

        query_plan_correctness_check(query_graph, vertex_orders[i]);

        std::vector<uint32_t> bn_cnt_list;
        uint32_t utility_value = query_plan_utility_value(query_graph, vertex_orders[i], bn_cnt_list);
        printf("query plan %d: utility value %d\n", i, utility_value);
        printf("\tvertex order: ");
        for (auto u : vertex_orders[i]) {
            printf("%3d ", u);
        }
        printf("\n");
        printf("\tbn:           ");
        for (auto cnt : bn_cnt_list) {
            printf("%3d ", cnt);
        }
        printf("\n");
        if (!node_orders.empty()) {
            printf("\tnode order: ");
            for (auto u : node_orders[i]) {
                printf("%3d ", u);
            }
            printf("\n");
        }

        if (utility_value > max_utility_value) {
            selected_vertex_order.clear();
            selected_vertex_order.emplace_back(vertex_orders[i]);
            max_utility_value = utility_value;
            selected_vertex_order_id = i;
        }
    }

    printf("selected vertex order id %d\n", selected_vertex_order_id);
    vertex_orders.swap(selected_vertex_order);
}

void query_plan_generator::query_plan_correctness_check(Graph *query_graph, std::vector<uint32_t> &vertex_order) {
    uint32_t n = query_graph->getVerticesCount();
    std::vector<bool> visited(n, false);
    if (vertex_order.size() != n) {
        printf("The number of vertices in the query plan is not correct %zu : %d.\n", vertex_order.size(), n);
        exit(-1);
    }
    for (auto u : vertex_order) {
        if (u >= n) {
            printf("Invalid vertex id %d in the query plan.\n", u);
            exit(-1);
        }

        if (visited[u]) {
            printf("Duplicate vertex id %d in the query plan.\n", u);
            exit(-1);
        }
        visited[u] = true;

        if (u != vertex_order[0]) {
            bool has_bn = false;
            uint32_t nbrs_cnt;
            const uint32_t *nbrs = query_graph->getVertexNeighbors(u, nbrs_cnt);

            for (uint32_t i = 0; i < nbrs_cnt; ++i) {
                if (visited[nbrs[i]]) {
                    has_bn = true;
                    break;
                }
            }

            if (!has_bn) {
                printf("Vertex id %d has no backward neighbors.\n", u);
                exit(-1);
            }
        }
    }
}

uint32_t query_plan_generator::query_plan_utility_value(Graph *query_graph, std::vector<uint32_t> &vertex_order,
                                                       std::vector<uint32_t> &bn_cn_list) {
    uint32_t n = query_graph->getVerticesCount();
    std::vector<bool> visited(n, false);

    for (auto u : vertex_order) {
        visited[u] = true;
        uint32_t nbrs_cnt;
        const uint32_t *nbrs = query_graph->getVertexNeighbors(u, nbrs_cnt);

        uint32_t bn_cnt = 0;
        for (uint32_t i = 0; i < nbrs_cnt; ++i) {
            if (visited[nbrs[i]]) {
                bn_cnt += 1;
            }
        }

        bn_cn_list.push_back(bn_cnt);
    }

    uint32_t sum = 0;

    for (uint32_t i = 0; i < n; ++i) {
        for (uint32_t j = 0; j <= i; ++j) {
            sum += bn_cn_list[j];
        }
    }

    return sum;
}

void query_plan_generator::update_extendable_vertex(Graph *query_graph, uint32_t u,
                                                    std::unordered_set<uint32_t> &extendable_vertex,
                                                    vector<bool> &visited_vertex) {
    uint32_t nbrs_cnt;
    const uint32_t* nbrs = query_graph->getVertexNeighbors(u, nbrs_cnt);

    for (uint32_t i = 0; i < nbrs_cnt; ++i) {
        uint32_t v = nbrs[i];
        if (!visited_vertex[v]) {
            extendable_vertex.insert(v);
        }
    }

    if (extendable_vertex.count(u) != 0) {
        extendable_vertex.erase(u);
    }
}

void query_plan_generator::greedy_expand(Graph *query_graph, catalog *storage, std::vector<uint32_t> &vertex_order,
                                         std::vector<bool> &visited_vertex,
                                         std::unordered_set<uint32_t> &extendable_vertex,
                                         uint32_t bn_cnt_threshold) {
    bool updated = true;

    while (updated) {
        updated = false;

        uint32_t max_bn_cnt = 0;
        uint32_t max_degree = 0;
        uint32_t min_candidate_vertex_num = std::numeric_limits<uint32_t>::max();
        uint32_t selected_vertex = 0;

        for (auto u : extendable_vertex) {
            uint32_t nbrs_cnt;
            const uint32_t* nbrs = query_graph->getVertexNeighbors(u, nbrs_cnt);

            uint32_t bn_cnt = 0;
            uint32_t local_candidate_vertex_num = storage->num_candidates_[u];
            uint32_t u_deg = query_graph->getVertexDegree(u);

            for (uint32_t i = 0; i < nbrs_cnt; ++i) {
                uint32_t v = nbrs[i];
                if (visited_vertex[v]) {
                    bn_cnt += 1;
                }
            }

            if (bn_cnt > max_bn_cnt || (bn_cnt == max_bn_cnt && local_candidate_vertex_num < min_candidate_vertex_num)
                || (bn_cnt == max_bn_cnt && local_candidate_vertex_num == min_candidate_vertex_num && u_deg > max_degree)) {
                max_bn_cnt = bn_cnt;
                max_degree = u_deg;
                min_candidate_vertex_num = local_candidate_vertex_num;
                selected_vertex = u;
            }
        }

        if (max_bn_cnt >= bn_cnt_threshold) {
            vertex_order.push_back(selected_vertex);
            visited_vertex[selected_vertex] = true;
            update_extendable_vertex(query_graph, selected_vertex, extendable_vertex, visited_vertex);
            updated = true;
        }
    }
}

void query_plan_generator::print_metrics() {
    printf("Ordering time (seconds): %.6f\n", NANOSECTOSEC(ordering_time_));
    printf("Nucleus decomposition time (seconds): %.6f\n", NANOSECTOSEC(nd_time_));
    printf("Traversal time (seconds): %.6f\n", NANOSECTOSEC(traversal_time_));
}