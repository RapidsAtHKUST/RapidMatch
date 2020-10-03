#include "scan.h"
#include <cstring>
#include <cereal/types/vector.hpp>
#include <vector>

void scan::execute(uint32_t src_label, uint32_t dst_label, edge_relation *relation, bool indexed) {
    if (indexed) {
        execute_with_index(src_label, dst_label, relation);
    }
    else {
        execute_without_index(src_label, dst_label, relation);
    }
}

void scan::execute(std::vector<uint32_t> &src_candidate_set, std::vector<uint32_t> &dst_candidate_set,
                   edge_relation *relation) {
    for (auto u : dst_candidate_set) {
        flag_[u] = true;
    }

    uint32_t edge_count = 0;

    for (auto u : src_candidate_set) {
        uint32_t dst_count;
        const uint32_t* dst_list = data_graph_->getVertexNeighbors(u, dst_count);

        for (uint32_t j = 0; j < dst_count; ++j) {
            uint32_t v = dst_list[j];

            if (flag_[v]) {
                buffer_[edge_count].vertices_[0] = u;
                buffer_[edge_count].vertices_[1] = v;
                edge_count += 1;
            }
        }
    }

    relation->size_ = edge_count;
    relation->edges_ = new edge[edge_count];
    memcpy(relation->edges_, buffer_, sizeof(edge) * edge_count);

    for (auto u : dst_candidate_set) {
        flag_[u] = false;
    }
}

void scan::execute_with_index(uint32_t src_label, uint32_t dst_label, edge_relation *relation) {
    uint64_t key = (uint64_t) src_label << 32 | dst_label;

    if (data_graph_->getEdgeIndex()->contains(key)) {
        auto edges = data_graph_->getEdgeIndex()->at(key);
        relation->size_ = edges->size();
        relation->edges_ = new edge[relation->size_];
        memcpy(relation->edges_, edges->data(), sizeof(edge) * relation->size_);
    }
    else {
        relation->size_ = 0;
    }
}

void scan::execute_without_index(uint32_t src_label, uint32_t dst_label, edge_relation *relation) {
    uint32_t edge_count = 0;

    uint32_t src_count;
    const uint32_t* src_list = data_graph_->getVerticesByLabel(src_label, src_count);

    for (uint32_t i = 0; i < src_count; ++i) {
        uint32_t src = src_list[i];

        uint32_t dst_count;
        const uint32_t* dst_list = data_graph_->getVertexNeighbors(src, dst_count);

        for (uint32_t j = 0; j < dst_count; ++j) {
            uint32_t dst = dst_list[j];

            if (data_graph_->getVertexLabel(dst) != dst_label)
                continue;

            buffer_[edge_count].vertices_[0] = src;
            buffer_[edge_count].vertices_[1] = dst;
            edge_count += 1;
        }
    }

    relation->size_ = edge_count;
    relation->edges_ = new edge[edge_count];
    memcpy(relation->edges_, buffer_, sizeof(edge) * edge_count);
}
