#ifndef SUBGRAPHMATCHING_FLAT_RELATION_H
#define SUBGRAPHMATCHING_FLAT_RELATION_H

#pragma once

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cassert>

/**
 * The default size limit for a relation is 4 GB.
 */
#define DEFAULT_FLAT_RELATION_CAPACITY (1024 * 1024 * 1024)

class flat_relation {
private:
    uint64_t capacity_;
    uint64_t size_;
    uint32_t length_;
    uint32_t * data_;

public:
    flat_relation(uint32_t length, uint64_t capacity = DEFAULT_FLAT_RELATION_CAPACITY)
            : capacity_(capacity), size_(0), length_(length) {
        data_ = (uint32_t*)malloc(capacity_ * sizeof(uint32_t));
    }

    ~flat_relation() {
        free(data_);
    }

    uint32_t get_tuple_length() {
        return length_;
    }

    uint64_t get_size() {
        return size_;
    }

    uint32_t* get_tuple(uint64_t offset) {
        assert(offset < size_);
        return data_ + (length_ * offset);
    }

    uint64_t push_back(uint32_t* tuple) {
        assert(length_ * (size_ + 1) <= capacity_);
        memcpy(data_ + length_ * size_, tuple, length_ * sizeof(uint32_t));
        size_ += 1;
        return size_ - 1;
    }

    uint64_t push_back(uint32_t* buffer, uint64_t count) {
        assert(size_ + count <= capacity_);
        memcpy(data_ + length_ * size_, buffer, count * sizeof(uint32_t));
        size_ += count;
        return size_ - 1;
    }

    uint32_t* push_back() {
        assert(length_ * (size_ + 1) <= capacity_);
        uint32_t* tuple = data_ + (length_ * size_);
        size_ += 1;
        return tuple;
    }

    uint64_t pop_pack() {
        assert(size_ > 0);
        size_ -= 1;
        return size_;
    }
};


#endif //SUBGRAPHMATCHING_FLAT_RELATION_H
