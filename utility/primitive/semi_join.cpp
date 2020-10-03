#include "semi_join.h"
#include <cstring>

uint32_t semi_join::execute(edge_relation *left, const uint32_t lkp, edge_relation *right, const uint32_t rkp) {
    // Initialize index.
    memset(index_, 0, sizeof(bool) * index_size_);

    uint32_t distinct_vertex_cnt = 0;
    // Build index.
    for (uint32_t i = 0; i < right->size_; ++i) {
        uint32_t key = right->edges_[i].vertices_[rkp];
        if (!index_[key]) {
            index_[key] = true;
            distinct_vertex_cnt += 1;
        }
    }

    // Semi join.
    uint32_t valid_edge_count = 0;
    for (uint32_t i = 0; i < left->size_; ++i) {
        uint32_t key = left->edges_[i].vertices_[lkp];

        if (index_[key]) {
            if (valid_edge_count != i) {
                left->edges_[valid_edge_count] = left->edges_[i];
            }
            valid_edge_count += 1;
        }
    }

    left->size_ = valid_edge_count;

    return distinct_vertex_cnt;
}
