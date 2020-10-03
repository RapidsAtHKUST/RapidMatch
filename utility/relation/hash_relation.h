#ifndef SUBGRAPHMATCHING_HASH_RELATION_H
#define SUBGRAPHMATCHING_HASH_RELATION_H

#include "../sparsepp/spp.h"
using spp::sparse_hash_map;

class hash_relation {
public:
    // Key: vertex id; Value: first is the number of children, and second is the offset.
    sparse_hash_map<uint32_t, std::pair<uint32_t, uint32_t>>* trie_;
    uint32_t cardinality_;
    uint32_t max_degree_;
    uint32_t* children_;
public:
    hash_relation() {
        trie_ = new sparse_hash_map<uint32_t, std::pair<uint32_t, uint32_t>>();
        cardinality_ = 0;
        max_degree_ = 0;
        children_ = nullptr;
    }

    ~hash_relation() {
        delete trie_;
        delete[] children_;
    }

    uint32_t get_children_num(uint32_t key) {
        auto iter = trie_->find(key);
        if (iter != trie_->end())
            return iter->second.first;
        return 0;
    }

    bool contains(uint32_t key) {
        return trie_->contains(key);
    }

    uint32_t* get_children(uint32_t key, uint32_t& count) {
        count = 0;
        auto iter = trie_->find(key);
        if (iter != trie_->end()) {
            count = iter->second.first;
            return children_ + iter->second.second;
        }
        return nullptr;
    }

    uint32_t get_size() {
        return static_cast<uint32_t>(trie_->size());
    }
    uint32_t get_cardinality() {
        return cardinality_;
    }

    uint64_t memory_cost() {
        if (get_size() == 0)
            return 0;

        uint64_t memory_cost = 0;
        uint64_t per_element_size = sizeof(uint32_t);
        memory_cost += get_cardinality() * per_element_size +  get_size() * 3 * per_element_size * trie_->load_factor();
        return memory_cost;
    }

    uint32_t max_degree() {
        return max_degree_;
    }
};

#endif //SUBGRAPHMATCHING_HASH_RELATION_H
