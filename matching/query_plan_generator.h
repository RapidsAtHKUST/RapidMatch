#ifndef SUBGRAPHMATCHING_QUERY_PLAN_GENERATOR_H
#define SUBGRAPHMATCHING_QUERY_PLAN_GENERATOR_H

#include "graph/graph.h"
#include "utility/relation/catalog.h"
#include "nucleus_decomposition/nd_interface.h"
#include <vector>
#include <unordered_set>

class query_plan_generator {
public:
    static double ordering_time_;
    static double nd_time_;
    static double traversal_time_;

private:
    static void construct_density_tree(std::vector<nd_tree_node>& density_tree, std::vector<nd_tree_node>& k12_tree,
                                       std::vector<nd_tree_node>& k23_tree, std::vector<nd_tree_node>& k34_tree);
    static void eliminate_node(std::vector<nd_tree_node>& density_tree, std::vector<nd_tree_node>& src_tree,
                               std::vector<nd_tree_node>& target_tree);
    static void merge_tree(std::vector<nd_tree_node>& density_tree);

    static void traversal_density_tree(Graph *query_graph, catalog *storage,
                                       std::vector<nd_tree_node> &density_tree,
                                       std::vector<std::vector<uint32_t>> &vertex_orders,
                                       std::vector<std::vector<uint32_t>> &node_orders);

    static void traversal_node(Graph *query_graph, catalog *storage, std::vector<nd_tree_node> &density_tree,
                                  std::vector<uint32_t> &vertex_order, std::vector<uint32_t> &node_order,
                                  vector<bool> &visited_vertex, std::vector<bool> &visited_node, nd_tree_node &cur_node,
                                  std::unordered_set<uint32_t> &extendable_vertex);

    static void greedy_expand(Graph *query_graph, catalog *storage, std::vector<uint32_t> &vertex_order,
                              std::vector<bool> &visited_vertex, std::unordered_set<uint32_t> &extendable_vertex,
                              uint32_t bn_cnt_threshold);

    static void update_extendable_vertex(Graph *query_graph, uint32_t u,
                                         std::unordered_set<uint32_t> &extendable_vertex,
                                         vector<bool> &visited_vertex);

    static double connectivity_common_neighbors(std::vector<bool> &visited_vertex, nd_tree_node &child_node);

    static void connectivity_shortest_path(Graph *query_graph, catalog *storage, std::vector<uint32_t> &vertex_order,
                                           std::vector<bool> &visited_vertex, nd_tree_node &cur_node,
                                           std::vector<uint32_t> &prev, std::vector<double> &dist);

    static void query_plan_correctness_check(Graph *query_graph, std::vector<uint32_t> &vertex_order);

    static uint32_t query_plan_utility_value(Graph *query_graph, std::vector<uint32_t> &vertex_order,
                                                std::vector<uint32_t> &bn_cn_list);

    static void print_density_tree(std::vector<nd_tree_node>& density_tree);

    public:
    static void generate_query_plan_for_test(Graph *query_graph, std::vector<uint32_t> &order);
    static void generate_query_plan_with_nd(Graph *query_graph, catalog *storage, std::vector<std::vector<uint32_t>>& vertex_orders);
    static void print_vertex_orders(Graph* query_graph, std::vector<std::vector<uint32_t>>& vertex_orders,
                                    std::vector<std::vector<uint32_t>>& node_orders);
    static void print_metrics();
};


#endif //SUBGRAPHMATCHING_QUERY_PLAN_GENERATOR_H
