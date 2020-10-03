#ifndef SUBGRAPHMATCHING_EXECUTION_TREE_NODE_H
#define SUBGRAPHMATCHING_EXECUTION_TREE_NODE_H

#include <vector>
#include <cstdint>
#include <sstream>
#include <chrono>
#include "utility/relation/flat_relation.h"
#include "utility/relation/catalog.h"
#include "utility/leapfrogtriejoin/leapfrogtriejoin.h"
#include "utility/hashjoin/nop_join.h"
#include "utility/cereal/archives/json.hpp"
#include "utility/cereal/types/vector.hpp"


#define OUTPUT_RELATION_ID(id) ((id) * 3 + 2)

enum join_operator_type {
    no_partition_hash_join = 0,
    leap_frog_trie_join = 1
};

class execution_tree_node {
protected:
    int id_;
    int parent_;
    join_operator_type operator_type_;
    std::vector<uint32_t> vertices_;

public:
    // Metrics
    uint64_t output_count_;
    double execution_time_;
    // double selectivity_;
    std::string cost_metrics_;
    std::string invalid_pr_count_;
    std::string fail_count_;

protected:
    virtual void initialize(catalog* catalog, std::vector<flat_relation *> &relations, uint64_t output_count_limit) = 0;
    virtual void execute_operator(std::vector<flat_relation *> &relations) = 0;
    virtual void clear() = 0;
    virtual std::string print_cost_metrics() = 0;
    virtual std::string print_invalid_pr_count() = 0;
    virtual std::string print_fail_count() = 0;
public:
    execution_tree_node(int id, int parent, join_operator_type type, std::vector<uint32_t>& vertices)
            : id_(id), parent_(parent), operator_type_(type), vertices_(vertices), output_count_(0), execution_time_(0) {}

    execution_tree_node() : output_count_(0), execution_time_(0) {}
    virtual ~execution_tree_node() {};

    join_operator_type get_operator_type() {
        return operator_type_;
    }

    int get_parent() {
        return parent_;
    }

    int get_id() {
        return id_;
    }

    std::vector<uint32_t>& get_vertices() {
        return vertices_;
    }

    bool is_root() {
        return parent_ == -1;
    }

    virtual bool is_leaf() = 0;

    uint64_t execute(catalog& catalog, std::vector<flat_relation*>& relations, uint64_t output_count_limit) {
        initialize(&catalog, relations, output_count_limit);

        auto start = std::chrono::high_resolution_clock::now();

        execute_operator(relations);

        auto end = std::chrono::high_resolution_clock::now();

        execution_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        cost_metrics_ = print_cost_metrics();
        invalid_pr_count_ = print_invalid_pr_count();
        fail_count_ = print_fail_count();

        clear();
        return output_count_;
    }
};

class no_partition_hash_join_node : public execution_tree_node {
protected:
    double hash_table_size_ratio_;
    int left_child_;
    int right_child_;
    flat_relation* left_relation_;
    flat_relation* right_relation_;
    flat_relation* output_;
    std::vector<uint32_t> left_child_vertices_;
    std::vector<uint32_t> right_child_vertices_;
    std::vector<uint32_t> left_child_key_position_;
    std::vector<uint32_t> right_child_key_position_;
    nop_join* nop_join_;

public:
    double left_child_cardinality_;
    double right_child_cardinality_;

protected:

public:
    no_partition_hash_join_node(int id, int parent, join_operator_type type, std::vector<uint32_t>& vertices, double hash_table_size_ratio,
                           int left_child, std::vector<uint32_t>& left_child_vertices, std::vector<uint32_t>& left_child_key_position,
                           int right_child, std::vector<uint32_t>& right_child_vertices, std::vector<uint32_t>& right_child_key_position)
            : execution_tree_node(id, parent, type, vertices), hash_table_size_ratio_(hash_table_size_ratio),
              left_child_(left_child), right_child_(right_child), left_relation_(nullptr), right_relation_(nullptr),
              output_(nullptr), left_child_vertices_(left_child_vertices), right_child_vertices_(right_child_vertices),
              left_child_key_position_(left_child_key_position), right_child_key_position_(right_child_key_position),
              nop_join_(nullptr)  {}

    no_partition_hash_join_node() : left_relation_(nullptr), right_relation_(nullptr), output_(nullptr), nop_join_(nullptr) {}
    double get_hash_table_size_ratio() {
        return hash_table_size_ratio_;
    }

    int get_left_child() {
        return left_child_;
    }

    std::vector<uint32_t>& get_left_child_vertices() {
        return left_child_vertices_;
    }

    std::vector<uint32_t>& get_left_child_key_position() {
        return left_child_key_position_;
    }

    int get_right_child() {
        return right_child_;
    }

    std::vector<uint32_t>& get_right_child_vertices() {
        return right_child_vertices_;
    }

    std::vector<uint32_t>& get_right_child_key_position() {
        return right_child_key_position_;
    }

    bool is_leaf () {
        return left_child_ == -1 && right_child_ == -1;
    }

    template <class Archive> void serialize(Archive& ar) {
        ar(cereal::make_nvp("Operator type", operator_type_),
           cereal::make_nvp("Hash table size ratio", hash_table_size_ratio_),
           cereal::make_nvp("Id", id_),
           cereal::make_nvp("Parent", parent_),
           cereal::make_nvp("Vertices", vertices_),
           cereal::make_nvp("Left child", left_child_),
           cereal::make_nvp("Left child vertices", left_child_vertices_),
           cereal::make_nvp("Left child key position", left_child_key_position_),
           cereal::make_nvp("Right child", right_child_),
           cereal::make_nvp("Right child vertices", right_child_vertices_),
           cereal::make_nvp("Right child key position", right_child_key_position_)
           );
    }
protected:
    void initialize(catalog *catalog, std::vector<flat_relation *> &relations, uint64_t output_count_limit) override;

    void execute_operator(std::vector<flat_relation *> &relations) override;

    void clear() override;

    std::string print_cost_metrics() override;

    std::string print_invalid_pr_count() override;

    std::string print_fail_count() override;
};

class leap_frog_trie_join_node : public execution_tree_node {
protected:
    int child_;
    std::vector<uint32_t> child_vertices_;
    flat_relation* input_relation_;
    flat_relation* output_relation_;
    leapfrogtriejoin* lftj_;
protected:

public:
    leap_frog_trie_join_node(int id, int parent, join_operator_type type, std::vector<uint32_t>& vertices,
            int child, std::vector<uint32_t>& child_vertices)
            : execution_tree_node(id, parent, type, vertices), child_(child), child_vertices_(child_vertices),
                input_relation_(nullptr), output_relation_(nullptr), lftj_(nullptr) {}

    leap_frog_trie_join_node() : input_relation_(nullptr), output_relation_(nullptr), lftj_(nullptr) {}

    int get_child() {
        return child_;
    }

    std::vector<uint32_t>& get_child_vertices() {
        return child_vertices_;
    }

    bool is_leaf() {
        return child_ == -1;
    }



    template <class Archive> void serialize(Archive& ar) {
        ar(cereal::make_nvp("Operator type", operator_type_),
           cereal::make_nvp("Id", id_),
           cereal::make_nvp("Parent", parent_),
           cereal::make_nvp("Vertices", vertices_),
           cereal::make_nvp("Child", child_),
           cereal::make_nvp("Child vertices", child_vertices_)
        );
    }
protected:
    void initialize(catalog *catalog, std::vector<flat_relation *> &relations, uint64_t output_count_limit) override;

    void execute_operator(std::vector<flat_relation *> &relations) override;

    void clear() override;

    std::string print_cost_metrics() override;

    std::string print_invalid_pr_count() override;

    std::string print_fail_count() override;

};

#endif //SUBGRAPHMATCHING_EXECUTION_TREE_NODE_H
