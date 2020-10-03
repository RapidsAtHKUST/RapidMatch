#include "execution_tree_node.h"
#include "utility/leapfrogtriejoin/leapfrogtriejoin.h"

void no_partition_hash_join_node::initialize(catalog *catalog, std::vector<flat_relation *> &relations,
                                             uint64_t output_count_limit) {
    assert(!is_leaf());
    left_relation_ = relations[OUTPUT_RELATION_ID(get_left_child())];
    right_relation_ = relations[OUTPUT_RELATION_ID(get_right_child())];

    // Set the left relation with smaller number of tuples.
//    if (left_relation_->get_size() > right_relation_->get_size()) {
//        std::swap(left_relation_, right_relation_);
//        left_child_vertices_.swap(right_child_vertices_);
//        left_child_key_position_.swap(right_child_key_position_);
//    }

    if (!is_root()) {
        output_ = new flat_relation(static_cast<uint32_t>(vertices_.size()));
    }

    if (left_child_key_position_.size() == 1) {
        nop_join_ = new nop_join(left_relation_, left_child_key_position_[0], right_relation_, right_child_key_position_[0],
                output_, HASH_TABLE_RATIO);
    }
    else if (left_child_key_position_.size() == 2) {
        nop_join_ = new nop_join(left_relation_, left_child_key_position_[0], left_child_key_position_[1], right_relation_,
                right_child_key_position_[0], right_child_key_position_[1], output_, HASH_TABLE_RATIO);
    }
    else {
        std::cerr << "Currently, we cannot support the pair-wise join with more than two vertices." << std::endl;
        exit(-1);
    }
}

void no_partition_hash_join_node::execute_operator(std::vector<flat_relation *> &relations) {
    output_count_ = nop_join_->execute();
    relations[OUTPUT_RELATION_ID(id_)] = output_;
}

void no_partition_hash_join_node::clear() {
    delete left_relation_;
    delete right_relation_;
    delete nop_join_;
}

std::string no_partition_hash_join_node::print_cost_metrics() {
    return "";
}

std::string no_partition_hash_join_node::print_invalid_pr_count() {
    return "";
}

std::string no_partition_hash_join_node::print_fail_count() {
    return "";
}

void leap_frog_trie_join_node::initialize(catalog *catalog, std::vector<flat_relation *> &relations, uint64_t output_count_limit) {
    if (!is_leaf()) {
        input_relation_ = relations[OUTPUT_RELATION_ID(child_)];
    }
    else {
#if RELATION_STRUCTURE == 0
        uint32_t u = vertices_[0];
        uint64_t count = catalog->get_num_candidates(u);
        input_relation_ = new flat_relation(1, count);
        input_relation_->push_back(catalog->get_candidates(u), count);
#elif RELATION_STRUCTURE == 1
        uint32_t u =vertices_[0];
        uint32_t v = vertices_[1];
        uint64_t count = catalog->hash_relations_[u][v].get_size();
        input_relation_ = new flat_relation(1, count);
        std::vector<uint32_t> buffer;
        buffer.reserve(count);

        for (auto& e : *catalog->hash_relations_[u][v].trie_)
            buffer.push_back(e.first);
        // Add this line for a fair comparison.
        std::sort(buffer.begin(), buffer.end());
        input_relation_->push_back(buffer.data(), count);
#elif RELATION_STRUCTURE == 2
        uint32_t u =vertices_[0];
        uint32_t v = vertices_[1];
        uint64_t count = catalog->trie_relations_[u][v].size_;
        input_relation_ = new flat_relation(1, count);
        input_relation_->push_back(catalog->trie_relations_[u][v].parent_, count);
#endif
    }

    // Determine output relation.
    if (!is_root()) {
        output_relation_ = new flat_relation(static_cast<uint32_t >(vertices_.size()));
    }

    lftj_ = new leapfrogtriejoin(input_relation_, output_relation_, vertices_.data(), static_cast<uint32_t >(vertices_.size()),
                                     catalog->query_graph_->get2CoreSize(), catalog, output_count_limit);
}

void leap_frog_trie_join_node::execute_operator(std::vector<flat_relation *> &relations) {
    output_count_ = lftj_->execute();
    relations[OUTPUT_RELATION_ID(id_)] = output_relation_;
}

void leap_frog_trie_join_node::clear() {
    delete input_relation_;
    delete lftj_;
}

std::string leap_frog_trie_join_node::print_cost_metrics() {
    std::stringstream ss;
    uint64_t total_pr_count = 0;
    uint64_t total_si_cost = 0;
    uint64_t total_si_count = 0;
    for (uint32_t i = static_cast<uint32_t>(child_vertices_.size()); i < vertices_.size(); ++i) {
        ss << "(" << vertices_[i] << "," << lftj_->intermediate_result_count_[i] << ","
           << lftj_->set_intersection_count_[i] << "," << lftj_->set_intersection_cost_[i] << ") ";
        total_pr_count += lftj_->intermediate_result_count_[i];
        total_si_cost += lftj_->set_intersection_cost_[i];
        total_si_count += lftj_->set_intersection_count_[i];
    }

    ss << "(Total," << total_pr_count << "," << total_si_count << "," << total_si_cost << "," << lftj_->fail_si_count_ << ")";
    return ss.str();
}

std::string leap_frog_trie_join_node::print_invalid_pr_count() {
    std::stringstream ss;
    ss << "(" << lftj_->invalid_core_pr_count_ + lftj_->invalid_leaf_pr_count_ << ","
        << lftj_->invalid_core_pr_count_ << "," <<lftj_->invalid_leaf_pr_count_ << ") ";
    return ss.str();
}

std::string leap_frog_trie_join_node::print_fail_count() {
    std::stringstream ss;
    ss << "(" << lftj_->iso_conflict_count_ + lftj_->fail_si_count_ << "," << lftj_->fail_si_count_ << ","
            << lftj_->iso_conflict_count_ << ") ";
    return ss.str();
}