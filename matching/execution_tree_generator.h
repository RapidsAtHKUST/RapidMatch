#ifndef SUBGRAPHMATCHING_EXECUTION_TREE_GENERATOR_H
#define SUBGRAPHMATCHING_EXECUTION_TREE_GENERATOR_H

#include "graph/graph.h"
#include "utility/relation/catalog.h"
#include "execution_tree.h"

class execution_tree_generator {
private:
    catalog* catalog_;

public:
    static void generate_wcoj_orders(const Graph* graph, std::vector<std::vector<uint32_t>>& orders, uint32_t num_limit);
    static execution_tree* generate_single_node_execution_tree(std::vector<uint32_t>& order);
};


#endif //SUBGRAPHMATCHING_EXECUTION_TREE_GENERATOR_H
