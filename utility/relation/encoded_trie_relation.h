#ifndef SUBGRAPHMATCHING_ENCODED_TRIE_RELATION_H
#define SUBGRAPHMATCHING_ENCODED_TRIE_RELATION_H

#include <cstdint>
#include <cassert>
#include <memory>
class encoded_trie_relation {
public:
    uint32_t size_;
    uint32_t max_degree_;
    uint32_t* offset_;
    uint32_t* children_;

public:
    encoded_trie_relation() : size_(0), max_degree_(0), offset_(nullptr), children_(nullptr) {}
    ~encoded_trie_relation() {
        delete[] offset_;
        delete[] children_;
    }

    uint32_t* get_children(uint32_t key, uint32_t& count) {
        assert(key < size_);
        count = offset_[key + 1] - offset_[key];
        return children_ + offset_[key];
    }

    uint32_t get_children_num(uint32_t key) {
        assert(key < size_);
        return offset_[key + 1] - offset_[key];
    }

    uint32_t get_cardinality() {
        return offset_[size_];
    }

    uint32_t get_size() {
        return size_;
    }

    uint64_t memory_cost() {
        if (size_ == 0)
            return 0;

        uint64_t memory_cost = 0;
        uint64_t per_element_size = sizeof(uint32_t);
        memory_cost += (size_ + 1) * per_element_size + get_cardinality() * per_element_size;
        return memory_cost;
    }

    uint32_t max_degree() {
        return max_degree_;
    }
};

#endif //SUBGRAPHMATCHING_ENCODED_TRIE_RELATION_H
