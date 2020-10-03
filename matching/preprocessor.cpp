#include "graphoperations.h"
#include "preprocessor.h"
#include "primitive/scan.h"
#include "primitive/nlf_filter.h"
#include "primitive/semi_join.h"
#include "primitive/projection.h"

void
preprocessor::execute(const Graph *query_graph, const Graph *data_graph, catalog *storage, bool enable_elimination) {
    initialize(query_graph, data_graph);
    scan_relation(storage);

    filter_time_ = 0;
#ifndef HOMOMORPHISM
    auto filter = new nlf_filter(query_graph, data_graph);
    auto start = std::chrono::high_resolution_clock::now();
    /**
     * NLF filter.
     */
    filter->execute(storage);
    auto end = std::chrono::high_resolution_clock::now();
    filter_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
#endif

    /**
     * -------------------------------------------------------
     * Collect statistics. The elapsed time is not counted.
     */
    for (uint32_t u = 0; u < vertices_count_; ++u) {
        for (ui v = u + 1; v < vertices_count_; ++v) {
            if (query_graph_->checkEdgeExistence(u, v)) {
                std::pair<uint32_t, uint32_t> key = std::make_pair(u, v);
                storage->catalog_info_[key].after_scan_ = storage->edge_relations_[u][v].size_;
            }
        }
    }

    /**
     * -------------------------------------------------------
     */

    if (enable_elimination) {
        generate_preprocess_plan();
        eliminate_dangling_tuples(storage);
    }

    /**
     * -------------------------------------------------------
     * Collect statistics. The elapsed time is not counted.
     */
    for (uint32_t u = 0; u < vertices_count_; ++u) {
        for (ui v = u + 1; v < vertices_count_; ++v) {
            if (query_graph_->checkEdgeExistence(u, v)) {
                std::pair<uint32_t, uint32_t> key = std::make_pair(u, v);
                storage->catalog_info_[key].after_eliminate_ = storage->edge_relations_[u][v].size_;
            }
        }
    }

    /**
     * -------------------------------------------------------
     */

    // TODO: Set max number vertices per vertex.
    for (uint32_t u = 0; u < vertices_count_; ++u) {
        for (uint32_t v = u + 1; v < vertices_count_; ++v) {
            if (query_graph->checkEdgeExistence(u, v)) {
                storage->max_num_candidates_per_vertex_ = std::max(storage->max_num_candidates_per_vertex_, storage->edge_relations_[u][v].size_);
            }
        }
    }

    preprocess_time_ = semi_join_time_ + scan_time_ + filter_time_;
}

void preprocessor::initialize(const Graph *query_graph, const Graph *data_graph) {
    query_graph_ = query_graph;
    data_graph_ = data_graph;
    vertices_count_ = query_graph_->getVerticesCount();
    non_core_vertices_count_ = vertices_count_ - query_graph_->get2CoreSize();
    degeneracy_ordering_ = new uint32_t[vertices_count_];
    vertices_index_ = new uint32_t[vertices_count_];

    if (non_core_vertices_count_ != 0) {
        non_core_vertices_parent_ = new uint32_t[non_core_vertices_count_];
        non_core_vertices_children_ = new uint32_t[non_core_vertices_count_];
        non_core_vertices_children_offset_ = new uint32_t[non_core_vertices_count_ + 1];
    }
    preprocess_time_ = 0;
    scan_time_ = 0;
    semi_join_time_ = 0;
    bottom_up_non_core_semi_join_time_ = 0;
    bottom_up_core_semi_join_time_ = 0;
    top_down_non_core_semi_join_time_ = 0;
    top_down_core_semi_join_time_ = 0;
}

void preprocessor::scan_relation(catalog *storage) {
    auto scan_operator = new scan(data_graph_);

    auto start = std::chrono::high_resolution_clock::now();
    for (uint32_t u = 0; u < vertices_count_; ++u) {
        uint32_t u_nbrs_cnt;
        const uint32_t* u_nbrs = query_graph_->getVertexNeighbors(u, u_nbrs_cnt);

        for (uint32_t i = 0; i < u_nbrs_cnt; ++i) {
            uint32_t v = u_nbrs[i];

            if (u > v)
                continue;

            uint32_t u_label = query_graph_->getVertexLabel(u);
            uint32_t v_label = query_graph_->getVertexLabel(v);
            scan_operator->execute(u_label, v_label, &(storage->edge_relations_[u][v]), true);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    scan_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    delete scan_operator;
}

void preprocessor::eliminate_dangling_tuples(catalog *storage) {
    auto semi_join_operator = new semi_join(data_graph_->getVerticesCount());
    auto start = std::chrono::high_resolution_clock::now();
    // Bottom-up semi-join along the degeneracy ordering for the non-core vertices.
    for (uint32_t i = 0; i < non_core_vertices_count_; ++i) {
        uint32_t u = degeneracy_ordering_[i];

        if (i == vertices_count_ - 1)
            break;

        uint32_t v = non_core_vertices_parent_[i];

        if (vertices_index_[v] < non_core_vertices_count_ && vertices_index_[v] < (vertices_count_ - 1)) {
            uint32_t w = non_core_vertices_parent_[vertices_index_[v]];

            // Left relation: R(v, w); Right relation: R(v, u); Semi join key: v.
            uint32_t lkp;
            edge_relation* l_relation = get_key_position_in_relation(v, w, storage, lkp);

            uint32_t rkp;
            edge_relation* r_relation = get_key_position_in_relation(v, u, storage, rkp);

            storage->num_candidates_[v] = semi_join_operator->execute(l_relation, lkp, r_relation, rkp);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    bottom_up_non_core_semi_join_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    // Perform the full reducer for the star rooted at each core vertex along the degeneracy ordering.
    for (uint32_t i = non_core_vertices_count_; i < vertices_count_; ++i) {
        uint32_t u = degeneracy_ordering_[i];
        uint32_t u_nbrs_cnt;
        const uint32_t* u_nbrs = query_graph_->getVertexNeighbors(u, u_nbrs_cnt);

        // Remove era along the order of neighbor sets.
        for (int j = 0; j < static_cast<int>(u_nbrs_cnt) - 1; ++j) {

            uint32_t v = u_nbrs[j];
            uint32_t w = u_nbrs[j + 1];

            // Left relation: R(u, w); Right relation: R(u, v); Semi join key u.
            uint32_t lkp;
            edge_relation* l_relation = get_key_position_in_relation(u, w, storage, lkp);
            uint32_t rkp;
            edge_relation* r_relation = get_key_position_in_relation(u, v, storage, rkp);

            storage->num_candidates_[u] = semi_join_operator->execute(l_relation, lkp, r_relation, rkp);
        }

        // Remove era along the reverse order of neighbor sets.
        for (int j = static_cast<int>(u_nbrs_cnt) - 1; j > 0; --j) {
            uint32_t v = u_nbrs[j - 1];
            uint32_t w = u_nbrs[j];

            // Left relation: R(u, v); Right relation: R(u, w); Semi join key u.
            uint32_t lkp;
            edge_relation* l_relation = get_key_position_in_relation(u, v, storage, lkp);
            uint32_t rkp;
            edge_relation* r_relation = get_key_position_in_relation(u, w, storage, rkp);

            storage->num_candidates_[u] = semi_join_operator->execute(l_relation, lkp, r_relation, rkp);
        }
    }

    end = std::chrono::high_resolution_clock::now();
    bottom_up_core_semi_join_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    // Perform the full reducer for the star rooted at each core vertex along the reverse of the degeneracy ordering.
    for (int i = static_cast<int>(vertices_count_) - 1; i >= static_cast<int>(non_core_vertices_count_); --i) {
        uint32_t u = degeneracy_ordering_[i];
        uint32_t u_nbrs_cnt;
        const uint32_t* u_nbrs = query_graph_->getVertexNeighbors(u, u_nbrs_cnt);

        // Remove era along the order of neighbor sets.
        for (int j = 0; j < static_cast<int>(u_nbrs_cnt) - 1; ++j) {

            uint32_t v = u_nbrs[j];
            uint32_t w = u_nbrs[j + 1];

            // Left relation: R(u, w); Right relation: R(u, v); Semi join key u.
            uint32_t lkp;
            edge_relation* l_relation = get_key_position_in_relation(u, w, storage, lkp);
            uint32_t rkp;
            edge_relation* r_relation = get_key_position_in_relation(u, v, storage, rkp);

            storage->num_candidates_[u] = semi_join_operator->execute(l_relation, lkp, r_relation, rkp);
        }

        // Remove era along the reverse order of neighbor sets.
        uint32_t cur_v = u_nbrs[u_nbrs_cnt - 1];
        uint32_t min_size = 0;
        if (u < cur_v) {
            min_size = storage->edge_relations_[u][cur_v].size_;
        }
        else {
            min_size = storage->edge_relations_[cur_v][u].size_;
        }

        for (int j = static_cast<int>(u_nbrs_cnt) - 1; j > 0; --j) {
            uint32_t v = u_nbrs[j - 1];
            uint32_t w = u_nbrs[j];

            // Left relation: R(u, v); Right relation: R(u, w); Semi join key u.
            uint32_t lkp;
            edge_relation* l_relation = get_key_position_in_relation(u, v, storage, lkp);
            uint32_t rkp;
            edge_relation* r_relation = get_key_position_in_relation(u, w, storage, rkp);

            storage->num_candidates_[u] = semi_join_operator->execute(l_relation, lkp, r_relation, rkp);

            if (l_relation->size_ < min_size) {
                min_size = l_relation->size_;
                cur_v = v;
            }
        }
    }

    end = std::chrono::high_resolution_clock::now();
    top_down_core_semi_join_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    // Top-down semi-join along the reverse of degeneracy ordering for the non-core vertices.
    for (int i = static_cast<int>(non_core_vertices_count_) - 1; i >= 0; --i) {
        // If all vertices are non-core vertices, then skip the last vertex in the degeneracy ordering.
        if (i == vertices_count_ - 1)
            continue;

        uint32_t u = degeneracy_ordering_[i];
        uint32_t v = non_core_vertices_parent_[i];

        uint32_t rkp;
        edge_relation* r_relation = get_key_position_in_relation(u, v, storage, rkp);

        uint32_t u_children_cnt = non_core_vertices_children_offset_[i + 1] - non_core_vertices_children_offset_[i];
        uint32_t* u_children = non_core_vertices_children_ + non_core_vertices_children_offset_[i];

        for (uint32_t j = 0; j < u_children_cnt; ++j) {
            uint32_t w = u_children[j];

            // Left relation: R(u, w); Right Relation: R(u, v); Semi join key: u.
            uint32_t lkp;
            edge_relation* l_relation = get_key_position_in_relation(u, w, storage, lkp);
            storage->num_candidates_[u] = semi_join_operator->execute(l_relation, lkp, r_relation, rkp);
        }
    }

    end = std::chrono::high_resolution_clock::now();
    top_down_non_core_semi_join_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    semi_join_time_ = bottom_up_non_core_semi_join_time_ + bottom_up_core_semi_join_time_ + top_down_non_core_semi_join_time_
            + top_down_core_semi_join_time_;
    delete semi_join_operator;
}

void preprocessor::generate_preprocess_plan() {
    GraphOperations::compute_degeneracy_order(query_graph_, degeneracy_ordering_);
    for (uint32_t i = 0; i < vertices_count_; ++i) {
        vertices_index_[degeneracy_ordering_[i]] = i;
    }

    if (non_core_vertices_count_ == 0)
        return;

    uint32_t children_offset = 0;
    for (uint32_t i = 0; i < non_core_vertices_count_; ++i) {
        uint32_t u = degeneracy_ordering_[i];
        uint32_t u_nbrs_cnt;
        const uint32_t* u_nbrs = query_graph_->getVertexNeighbors(u, u_nbrs_cnt);
        non_core_vertices_children_offset_[i] = children_offset;

        for (uint32_t j = 0; j < u_nbrs_cnt; ++j) {
            uint32_t v = u_nbrs[j];

            if (vertices_index_[v] < i) {
                non_core_vertices_children_[children_offset++] = v;
            }
            else {
                non_core_vertices_parent_[i] = v;
            }
        }
    }

    non_core_vertices_children_offset_[non_core_vertices_count_] = children_offset;
}

edge_relation *preprocessor::get_key_position_in_relation(uint32_t u, uint32_t v, catalog *storage, uint32_t &kp) {
    kp = 0;

    if (u > v) {
        kp = 1;
        std::swap(u, v);
    }
    return &(storage->edge_relations_[u][v]);
}

void preprocessor::print_metrics() {
    printf("Preprocessing time (seconds): %.6f\n", NANOSECTOSEC(semi_join_time_ + scan_time_ + filter_time_));
    printf("Filter time (seconds): %.6f\n", NANOSECTOSEC(filter_time_));
    printf("Scan time (seconds): %.6f\n", NANOSECTOSEC(scan_time_));
    printf("Semi-join time (seconds): %.6f\n", NANOSECTOSEC(semi_join_time_));
    printf("Bottom-up semi-join time on non-core relations (seconds): %.6f\n", NANOSECTOSEC(bottom_up_non_core_semi_join_time_));
    printf("Bottom-up semi-join time on core relations (seconds): %.6f\n", NANOSECTOSEC(bottom_up_core_semi_join_time_));
    printf("Top-down semi-join time on non-core relations (seconds): %.6f\n", NANOSECTOSEC(top_down_non_core_semi_join_time_));
    printf("Top-down semi-join time on core relations (seconds): %.6f\n", NANOSECTOSEC(top_down_core_semi_join_time_));
}
