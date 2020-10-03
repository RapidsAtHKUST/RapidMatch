#include "execution_tree_generator.h"

void execution_tree_generator::generate_wcoj_orders(const Graph *graph, std::vector<std::vector<uint32_t>> &orders,
                                                    uint32_t num_limit) {
    uint32_t n = graph->getVerticesCount();
    std::vector<bool> visited(n, false);
    std::vector<uint32_t> idx(n);
    std::vector<uint32_t> order;

    for (uint32_t u = 0; u < n; ++u) {
        order.push_back(u);
        visited[u] = true;

        uint32_t cur_depth = 1;
        idx[cur_depth] = 0;

        while (true) {
            while (idx[cur_depth] < n) {
                uint32_t v = idx[cur_depth]++;

                if (visited[v])
                    continue;

                bool valid = false;
                for (auto temp : order) {
                    if (graph->checkEdgeExistence(v, temp)) {
                        valid = true;
                        break;
                    }
                }

                if (!valid)
                    continue;

                order.push_back(v);
                visited[v] = true;

                if (order.size() == n) {
                    orders.emplace_back(order);
                    visited[v] = false;
                    order.pop_back();

                    if (orders.size() >= num_limit)
                        return;
                }
                else {
                    cur_depth += 1;
                    idx[cur_depth] = 0;
                }
            }

            cur_depth -= 1;
            visited[order.back()] = false;
            order.pop_back();
            if (cur_depth == 0)
                break;
        }
    }
}

execution_tree *execution_tree_generator::generate_single_node_execution_tree(std::vector<uint32_t> &order) {
    std::vector<uint32_t> child;
    std::vector<execution_tree_node*> nodes;
    auto node = new leap_frog_trie_join_node(0, -1, join_operator_type::leap_frog_trie_join, order, -1, child);
    nodes.push_back(node);
    return new execution_tree(nodes);
}
