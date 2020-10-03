#include "encoder.h"
#include "primitive/projection.h"

void encoder::convert_to_hash_relation(catalog *storage, uint32_t u, uint32_t v) {
    // We assume that the relation is ordered.
    uint32_t src = std::min(u, v);
    uint32_t dst = std::max(u, v);
    edge_relation& target_edge_relation = storage->edge_relations_[src][dst];
    hash_relation& target_hash_relation1 = storage->hash_relations_[u][v];

    edge* edges = target_edge_relation.edges_;
    uint32_t edge_size = target_edge_relation.size_;

    assert(edge_size > 0);

    uint32_t u_key = 0;
    uint32_t v_key = 1;

    if (src != u) {
        std::swap(u_key, v_key);
        // Sort the target edge relation.
        std::sort(edges, edges + edge_size, [](edge& l, edge& r)-> bool {
            if (l.vertices_[1] == r.vertices_[1])
                return l.vertices_[0] < r.vertices_[0];
            return l.vertices_[1] < r.vertices_[1];
        });
    }

    target_hash_relation1.children_ = new uint32_t[edge_size];

    uint32_t offset = 0;
    uint32_t local_degree = 0;
    uint32_t prev_u0 = max_vertex_id_ + 1;

    for (uint32_t i = 0; i < edge_size; ++i) {
        uint32_t u0 = edges[i].vertices_[u_key];
        uint32_t u1 = edges[i].vertices_[v_key];
        if (u0 != prev_u0 ) {
            if (prev_u0 != max_vertex_id_ + 1)
                target_hash_relation1.trie_->emplace(prev_u0, std::make_pair(local_degree, offset));

            offset += local_degree;

            if (local_degree > target_hash_relation1.max_degree_) {
                target_hash_relation1.max_degree_ = local_degree;
            }

            local_degree = 0;
            prev_u0 = u0;
        }

        target_hash_relation1.children_[offset + local_degree] = u1;
        local_degree += 1;
    }

    target_hash_relation1.cardinality_ = edge_size;
    target_hash_relation1.trie_->emplace(prev_u0, std::make_pair(local_degree, offset));
    if (local_degree > target_hash_relation1.max_degree_) {
        target_hash_relation1.max_degree_ = local_degree;
    }
}

void encoder::convert_to_trie_relation(catalog *storage, uint32_t u, uint32_t v) {
    // We assume that the relation is ordered.
    uint32_t src = std::min(u, v);
    uint32_t dst = std::max(u, v);

    edge_relation& target_edge_relation = storage->edge_relations_[src][dst];
    trie_relation& trie_relation1 = storage->trie_relations_[u][v];

    edge* edges = target_edge_relation.edges_;
    uint32_t edge_size = target_edge_relation.size_;

    assert(edge_size > 0);

    uint32_t u_key = 0;
    uint32_t v_key = 1;

    if (src != u) {
        std::swap(u_key, v_key);
        // Sort the target edge relation.
        std::sort(edges, edges + edge_size, [](edge& l, edge& r)-> bool {
            if (l.vertices_[1] == r.vertices_[1])
                return l.vertices_[0] < r.vertices_[0];
            return l.vertices_[1] < r.vertices_[1];
        });
    }

    trie_relation1.children_ = new uint32_t[edge_size];

    uint32_t offset = 0;
    uint32_t prev_u0 = max_vertex_id_ + 1;

    for (uint32_t i = 0; i < edge_size; ++i) {
        uint32_t u0 = edges[i].vertices_[u_key];
        uint32_t u1 = edges[i].vertices_[v_key];
        if (u0 != prev_u0 ) {
            temp_buffer1_[trie_relation1.size_] = u0;
            temp_buffer2_[trie_relation1.size_++] = offset;

            uint32_t local_degree = offset - temp_buffer2_[trie_relation1.size_ - 1];
            if (local_degree > trie_relation1.max_degree_) {
                trie_relation1.max_degree_ = local_degree;
            }

            prev_u0 = u0;
        }

        trie_relation1.children_[offset++] = u1;
    }

    temp_buffer2_[trie_relation1.size_] = offset;

    uint32_t local_degree = offset - temp_buffer2_[trie_relation1.size_ - 1];
    if (local_degree > trie_relation1.max_degree_) {
        trie_relation1.max_degree_ = local_degree;
    }

    trie_relation1.parent_ = new uint32_t[trie_relation1.size_];
    trie_relation1.offset_ = new uint32_t[trie_relation1.size_ + 1];

    memcpy(trie_relation1.parent_, temp_buffer1_, sizeof(uint32_t) * trie_relation1.size_);
    memcpy(trie_relation1.offset_, temp_buffer2_, sizeof(uint32_t) * (trie_relation1.size_ + 1));
}

void encoder::convert_to_encoded_relation(catalog *storage) {
    uint32_t core_vertices_cnt = query_graph_->get2CoreSize();

    auto projection_operator = new projection(max_vertex_id_);
    auto start = std::chrono::high_resolution_clock::now();
    for (uint32_t i = 0; i < core_vertices_cnt || i == 0; ++i) {
        uint32_t u = join_plan_[i];
        uint32_t nbr_cnt;
        const uint32_t* nbrs = query_graph_->getVertexNeighbors(u, nbr_cnt);
        for (uint32_t j = 0; j < nbr_cnt; ++j) {
            uint32_t v = nbrs[j];
            uint32_t src = u;
            uint32_t dst = v;
            uint32_t kp = 0;
            if (src > dst) {
                std::swap(src, dst);
                kp = 1;
            }

            projection_operator->execute(&storage->edge_relations_[src][dst], kp, storage->candidate_sets_[u], storage->num_candidates_[u]);
            break;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    projection_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    delete projection_operator;

    uint32_t n = query_graph_->getVerticesCount();
    for (uint32_t i = 1; i < n; ++i) {
        uint32_t u = join_plan_[i];
        for (uint32_t j = 0; j < i; ++j) {
            uint32_t bn = join_plan_[j];
            if (query_graph_->checkEdgeExistence(bn, u)) {
                if (i < core_vertices_cnt) {
                    convert_to_encoded_relation(storage, bn, u);
                }
                else {
                    convert_to_hash_relation(storage, bn, u);
                }
            }
        }
    }
}

void encoder::execute(catalog *storage, RelationStructure relation_type, bool enable_sparsebp, uint32_t *join_plan) {
    auto start = std::chrono::high_resolution_clock::now();
    join_plan_ = join_plan;
    switch (relation_type) {
        case RelationStructure::TrieRelation:
            convert_to_trie_relation(storage);
            break;
        case RelationStructure::HashRelation:
            convert_to_hash_relation(storage);
            break;
        case RelationStructure::EncodedTrieRelation:
            convert_to_encoded_relation(storage);
            break;
        default:
            std::cerr << "The relation structure is not supported." << std::endl;
            exit(-1);
    }

    auto end = std::chrono::high_resolution_clock::now();
    build_relation_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    if (enable_sparsebp) {
        start = std::chrono::high_resolution_clock::now();
        switch (relation_type) {
            case RelationStructure::TrieRelation:
                convert_trie_relation_to_sparse_bitmap(storage);
                break;
            case RelationStructure::HashRelation:
                convert_hash_relation_to_sparse_bitmap(storage);
                break;
            case RelationStructure::EncodedTrieRelation:
                convert_encoded_relation_to_sparse_bitmap(storage);
                break;
        }
        end = std::chrono::high_resolution_clock::now();
        convert_to_sparse_bp_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    }

    encoding_time_ = build_relation_time_ + convert_to_sparse_bp_time_;

    /**
     * -------------------------------------------------------
     * Collect statistics. The elapsed time is not counted.
     */
    for (uint32_t u = 0; u < query_graph_->getVerticesCount(); ++u) {
        for (ui v = u + 1; v < query_graph_->getVerticesCount(); ++v) {
            if (query_graph_->checkEdgeExistence(u, v)) {
                std::pair<uint32_t, uint32_t> key = std::make_pair(u, v);
                auto& meta_info = storage->catalog_info_[key];
                meta_info.max_degree_ = storage->get_max_degree(u, v);
                meta_info.reverse_max_degree_ = storage->get_max_degree(v, u);
            }
        }
    }
    /**
     * -------------------------------------------------------
     */
}

void encoder::convert_to_encoded_relation(catalog *storage, uint32_t u, uint32_t v) {
    uint32_t src = std::min(u, v);
    uint32_t dst = std::max(u, v);
    edge_relation& target_edge_relation = storage->edge_relations_[src][dst];
    edge* edges = target_edge_relation.edges_;
    uint32_t edge_size = target_edge_relation.size_;
    assert(edge_size > 0);

    uint32_t v_candidates_cnt = storage->get_num_candidates(v);
    uint32_t* v_candidates = storage->get_candidates(v);

    for (uint32_t i = 0; i < v_candidates_cnt; ++i) {
        uint32_t v_candidate = v_candidates[i];
        temp_buffer1_[v_candidate] = i + 1;
    }

    uint32_t u_p = 0;
    uint32_t v_p = 1;
    if (u > v) {
        // Sort R(v, u) by u.
        std::sort(edges, edges + edge_size, [](edge& l, edge& r) -> bool {
            if (l.vertices_[1] == r.vertices_[1])
                return l.vertices_[0] < r.vertices_[0];
            return l.vertices_[1] < r.vertices_[1];
        });
        u_p = 1;
        v_p = 0;
    }

    encoded_trie_relation& target_encoded_trie_relation = storage->encoded_trie_relations_[u][v];
    uint32_t edge_cnt = edge_size;
    uint32_t u_candidates_cnt = storage->get_num_candidates(u);
    uint32_t* u_candidates = storage->get_candidates(u);
    target_encoded_trie_relation.size_ = u_candidates_cnt;
    target_encoded_trie_relation.offset_ = new uint32_t[u_candidates_cnt + 1];
    target_encoded_trie_relation.children_ = new uint32_t[edge_size];

    uint32_t offset = 0;
    uint32_t edge_index = 0;

    for (uint32_t i = 0; i < u_candidates_cnt; ++i) {
        uint32_t u_candidate = u_candidates[i];
        target_encoded_trie_relation.offset_[i] = offset;
        uint32_t local_degree = 0;
        while (edge_index < edge_cnt) {
            uint32_t u0 = edges[edge_index].vertices_[u_p];
            uint32_t v0 = edges[edge_index].vertices_[v_p];
            if (u0 == u_candidate) {
                if (temp_buffer1_[v0] > 0) {
                    target_encoded_trie_relation.children_[offset + local_degree] = temp_buffer1_[v0] - 1;
                    local_degree += 1;
                }
            }
            else if (u0 > u_candidate) {
                break;
            }

            edge_index += 1;
        }

        offset += local_degree;

        if (local_degree > target_encoded_trie_relation.max_degree_) {
            target_encoded_trie_relation.max_degree_ = local_degree;
        }
    }

    target_encoded_trie_relation.offset_[u_candidates_cnt] = offset;

    for (uint32_t i = 0; i < v_candidates_cnt; ++i) {
        uint32_t v_candidate = v_candidates[i];
        temp_buffer1_[v_candidate] = 0;
    }
}

void encoder::convert_trie_relation_to_sparse_bitmap(catalog *storage) {
    uint32_t core_vertices_cnt = query_graph_->get2CoreSize();

    for (uint32_t i = 1; i < core_vertices_cnt; ++i) {
        uint32_t u = join_plan_[i];

        for (uint32_t j = 0; j < i; ++j) {
            uint32_t bn = join_plan_[j];
            if (query_graph_->checkEdgeExistence(u, bn)) {
                storage->bsr_relations_[bn][u].load(storage->trie_relations_[bn][u].get_size(),
                                                   storage->trie_relations_[bn][u].parent_,
                                                   storage->trie_relations_[bn][u].offset_,
                                                   storage->trie_relations_[bn][u].children_,
                                                   storage->max_num_candidates_per_vertex_);
            }
        }
    }
}

void encoder::convert_encoded_relation_to_sparse_bitmap(catalog *storage) {
    uint32_t core_vertices_cnt = query_graph_->get2CoreSize();

    for (uint32_t i = 1; i < core_vertices_cnt; ++i) {
        uint32_t u = join_plan_[i];

        for (uint32_t j = 0; j < i; ++j) {
            uint32_t bn = join_plan_[j];
            if (query_graph_->checkEdgeExistence(u, bn)) {
                storage->bsr_relations_[bn][u].load(storage->encoded_trie_relations_[bn][u].get_size(),
                                                   storage->encoded_trie_relations_[bn][u].offset_,
                                                   storage->encoded_trie_relations_[bn][u].offset_,
                                                   storage->encoded_trie_relations_[bn][u].children_,
                                                   storage->max_num_candidates_per_vertex_, true);
            }
        }
    }
}

void encoder::convert_hash_relation_to_sparse_bitmap(catalog *storage) {
    uint32_t core_vertices_cnt = query_graph_->get2CoreSize();

    for (uint32_t i = 1; i < core_vertices_cnt; ++i) {
        uint32_t u = join_plan_[i];

        for (uint32_t j = 0; j < i; ++j) {
            uint32_t bn = join_plan_[j];
            if (query_graph_->checkEdgeExistence(u, bn)) {
                storage->bsr_relations_[bn][u].load(storage->hash_relations_[bn][u].trie_,
                                                   storage->hash_relations_[bn][u].children_,
                                                   storage->max_num_candidates_per_vertex_);
            }
        }
    }
}

void encoder::convert_to_hash_relation(catalog *storage) {
    uint32_t n = query_graph_->getVerticesCount();

    for (uint32_t i = 1; i < n; ++i) {
        uint32_t u = join_plan_[i];

        for (uint32_t j = 0; j < i; ++j) {
            uint32_t bn = join_plan_[j];

            if (query_graph_->checkEdgeExistence(bn, u)) {
                convert_to_hash_relation(storage, bn, u);
            }
        }
    }
}

void encoder::convert_to_trie_relation(catalog *storage) {
    uint32_t n = query_graph_->getVerticesCount();

    for (uint32_t i = 1; i < n; ++i) {
        uint32_t u = join_plan_[i];

        for (uint32_t j = 0; j < i; ++j) {
            uint32_t bn = join_plan_[j];
            if (query_graph_->checkEdgeExistence(bn, u)) {
                convert_to_trie_relation(storage, bn, u);
            }
        }
    }
}

void encoder::print_metrics() {
    printf("Encoding time (seconds): %.6f\n", NANOSECTOSEC(encoding_time_));
    printf("Build relation time (seconds): %.6f\n", NANOSECTOSEC(build_relation_time_));
    printf("Project time (seconds): %.6f\n", NANOSECTOSEC(projection_time_));
    printf("Convert to sparsebp time (seconds): %.6f\n", NANOSECTOSEC(convert_to_sparse_bp_time_));
}


