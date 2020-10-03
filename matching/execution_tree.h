#ifndef SUBGRAPHMATCHING_EXECUTION_TREE_H
#define SUBGRAPHMATCHING_EXECUTION_TREE_H

#include <vector>
#include "execution_tree_node.h"
class execution_tree {
private:
    std::vector<execution_tree_node*> tree_nodes_;
    std::vector<flat_relation*> relations_;

private:
    void import_plan(const std::string plan);

public:
    execution_tree(std::vector<execution_tree_node*>& tree_nodes) : tree_nodes_(tree_nodes) {
        relations_.resize(tree_nodes_.size() * 3, nullptr);
    }

    execution_tree(const std::string& plan) {
        import_plan(plan);
        relations_.resize(tree_nodes_.size() * 3, nullptr);
    }

    std::string export_plan();

    ~execution_tree() {
        for (auto node : tree_nodes_) {
            delete node;
        }

        for (auto relation : relations_) {
            delete relation;
        }
    }

    uint64_t execute(catalog& catalog, uint64_t output_count_limit);

    double get_execution_time() {
        double execution_time = 0;
        for (auto& node : tree_nodes_) {
            execution_time += node->execution_time_;
        }
        return execution_time;
    }

    uint64_t get_output_count() {
        return tree_nodes_.back()->output_count_;
    }

    std::string get_cost_metrics() {
        std::string str;
        for (auto& node : tree_nodes_) {
            str += node->cost_metrics_;
        }
        return str;
    }

    std::string get_invalid_pr_count() {
        std::string str;
        for (auto& node : tree_nodes_) {
            str += node->invalid_pr_count_;
        }
        return str;
    }

    std::string get_fail_count() {
        std::string str;
        for (auto& node : tree_nodes_) {
            str += node->fail_count_;
        }
        return str;
    }
};


#endif //SUBGRAPHMATCHING_EXECUTION_TREE_H
