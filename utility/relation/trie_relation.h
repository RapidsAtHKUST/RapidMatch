#ifndef SUBGRAPHMATCHING_TRIE_RELATION_H
#define SUBGRAPHMATCHING_TRIE_RELATION_H

#include "../primitive/search.h"

class trie_relation {
public:
    uint32_t size_;
    uint32_t max_degree_;
    uint32_t* parent_;
    uint32_t* offset_;
    uint32_t* children_;

public:
    trie_relation() : size_(0), max_degree_(0), parent_(nullptr), offset_(nullptr), children_(nullptr) {}
    ~trie_relation() {
        delete[] parent_;
        delete[] offset_;
        delete[] children_;
    }

    uint32_t get_children_num(uint32_t key) {
        uint32_t index = search::binary_search(parent_, 0, size_, key);
        if (index == size_) {
            return 0;
        }

        return offset_[index + 1] - offset_[index];
    }

    uint32_t get_size() {
        return size_;
    }

    bool contains(uint32_t key) {
        uint32_t index = search::binary_search(parent_, 0, size_, key);
        return index != size_;
    }

    uint32_t get_cardinality() {
        return offset_[size_];
    }

    uint32_t* get_children(uint32_t key, uint32_t& count) {
        count = 0;
        uint32_t index = search::binary_search(parent_, 0, size_, key);
        if (index == size_) {
            return nullptr;
        }

        count = offset_[index + 1] - offset_[index];
        return children_ + offset_[index];
    }

    uint64_t memory_cost() {
        if (size_ == 0)
            return 0;

        uint64_t memory_cost = 0;
        uint64_t per_element_size = sizeof(uint32_t);
        memory_cost += get_cardinality() * per_element_size + get_size() * 2 * per_element_size;
        return memory_cost;
    }

    uint32_t max_degree() {
        return max_degree_;
    }
};


#endif //SUBGRAPHMATCHING_TRIE_RELATION_H
