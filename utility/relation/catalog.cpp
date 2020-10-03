#include "catalog.h"
#include "pretty_print.h"
#include <vector>

uint64_t catalog::memory_cost() {
    uint64_t memory_cost = 0;
    uint64_t per_element_size = sizeof(uint32_t);

#if RELATION_STRUCTURE == 0
    for (uint32_t i = 0; i < num_sets_; ++i) {
        memory_cost += num_candidates_[i] * per_element_size;
    }
#endif

    for (uint32_t i = 0; i < num_sets_; ++i) {
        for (uint32_t j = 0; j < num_sets_; ++j) {
            memory_cost += hash_relations_[i][j].memory_cost();
            memory_cost += trie_relations_[i][j].memory_cost();
            memory_cost += encoded_trie_relations_[i][j].memory_cost();
        }
    }

    return memory_cost;
}

void catalog::print_metrics() {
    uint32_t core_edges_cnt = 0;
    uint32_t tree_edges_cnt = 0;
    uint32_t leaf_edges_cnt = 0;
    for (auto& meta_info : catalog_info_) {
        if (meta_info.second.type == EdgeType::CoreEdge)
            core_edges_cnt += 1;
        else if (meta_info.second.type == EdgeType::TreeEdge)
            tree_edges_cnt += 1;
        else
            leaf_edges_cnt += 1;
    }

    printf("Index Info: CoreRelation(%u), TreeRelation(%u), LeafRelation(%u)\n", core_edges_cnt, tree_edges_cnt, leaf_edges_cnt);
    printf("Table Header: Relation Name, Relation Size after Scan, Relation Size after Elimination, Max Degree, Reverse Max Degree\n");
    for (auto& meta_info : catalog_info_) {
        std::string ss;
        if (meta_info.second.type == EdgeType::CoreEdge)
            ss = "CoreRelation";
        else if (meta_info.second.type == EdgeType::TreeEdge)
            ss = "TreeRelation";
        else
            ss = "LeafRelation";

        printf("%s %u-%u: %u, %u, %u, %u\n", ss.c_str(), meta_info.first.first, meta_info.first.second,
               meta_info.second.after_scan_, meta_info.second.after_eliminate_,
               meta_info.second.max_degree_, meta_info.second.reverse_max_degree_);
    }
}

void catalog::print_cardinality() {
    std::vector<std::pair<std::pair<uint32_t, uint32_t >, uint32_t>> core_edges;
    std::vector<std::pair<std::pair<uint32_t, uint32_t >, uint32_t>> tree_edges;
    std::vector<std::pair<std::pair<uint32_t, uint32_t >, uint32_t>> leaf_edges;

#ifdef SPARSE_BITMAP
    std::vector<std::pair<std::pair<uint32_t, uint32_t >, std::pair<double, uint32_t>>> core_edges_compact_ratio;
#endif
    double sum = 0;
    for (uint32_t i = 0; i < num_sets_; ++i) {
        uint32_t begin_vertex = i;
        for (ui j = i + 1; j < num_sets_; ++j) {
            uint32_t end_vertex = j;
            uint32_t edge_type = 0;
            if (query_graph_->checkEdgeExistence(begin_vertex, end_vertex)) {
                uint32_t cardinality = 0;

                if (hash_relations_[begin_vertex][end_vertex].memory_cost() != 0) {
                    cardinality = hash_relations_[begin_vertex][end_vertex].get_cardinality();
                }
                else if (encoded_trie_relations_[begin_vertex][end_vertex].memory_cost() != 0) {
                    cardinality = encoded_trie_relations_[begin_vertex][end_vertex].get_cardinality();
                }
                else if (trie_relations_[begin_vertex][end_vertex].memory_cost() != 0) {
                    cardinality = trie_relations_[begin_vertex][end_vertex].get_cardinality();
                }
                else {
                    std::cerr << "Relation has no data." << std::endl;
                    exit(-1);
                }

                sum += cardinality;
                if (query_graph_->getCoreValue(begin_vertex) > 1 && query_graph_->getCoreValue(end_vertex) > 1) {
                    core_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
                }
                else if (query_graph_->getVertexDegree(begin_vertex) == 1 || query_graph_->getVertexDegree(end_vertex) == 1) {
                    leaf_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
                    edge_type = 1;
                }
                else {
                    tree_edges.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), cardinality));
                    edge_type = 2;
                }

#ifdef SPARSE_BITMAP
                if (edge_type == 0) {
                    double ratio = 0;
                    uint32_t small_degree_count = 0;

#if RELATION_STRUCTURE == 0
                    uint32_t size = encoded_trie_relations_[begin_vertex][end_vertex].get_size();
                    for (uint32_t key = 0; key < size; ++key) {
                        uint32_t uncompact_size = encoded_trie_relations_[begin_vertex][end_vertex].get_children_num(key);
                        uint32_t compact_size = bsr_relations_[begin_vertex][end_vertex].bsrs[key].size_;
#elif RELATION_STRUCTURE == 1
                    uint32_t size = hash_relations_[begin_vertex][end_vertex].get_size();
                    for (auto& e : *hash_relations_[begin_vertex][end_vertex].trie_) {
                        uint32_t uncompact_size = e.second.first;
                        uint32_t compact_size = bsr_relations_[begin_vertex][end_vertex].bsrs[e.first].size_;
#elif RELATION_STRUCTURE == 2
                    uint32_t size = trie_relations_[begin_vertex][end_vertex].get_size();
                    for (uint32_t k = 0; k < size; ++k) {
                        uint32_t key = trie_relations_[begin_vertex][end_vertex].parent_[k];
                        uint32_t uncompact_size = trie_relations_[begin_vertex][end_vertex].get_children_num(key);
                        uint32_t compact_size = bsr_relations_[begin_vertex][end_vertex].bsrs[key].size_;
#endif
                        if (uncompact_size >= 4)
                            ratio += compact_size / (double) uncompact_size;
                        else
                            small_degree_count += 1;
                    }
                    ratio /= size;
                    core_edges_compact_ratio.emplace_back(std::make_pair(std::make_pair(begin_vertex, end_vertex), std::make_pair(ratio, small_degree_count)));
                }
#endif
            }
        }
    }
    printf("Max Degree: %u\n", max_num_candidates_per_vertex_);
    printf("Candidate Set Info:\n");
    for (uint32_t u = 0; u < num_sets_; ++u) {
        printf("CandidateSet %u: %u\n", u, num_candidates_[u]);
    }

    printf("Index Info: CoreTable(%zu), TreeTable(%zu), LeafTable(%zu)\n", core_edges.size(), tree_edges.size(), leaf_edges.size());

    for (auto table_info : core_edges) {
        printf("CoreTable %u-%u: %u\n", table_info.first.first, table_info.first.second, table_info.second);
    }

    for (auto table_info : tree_edges) {
        printf("TreeTable %u-%u: %u\n", table_info.first.first, table_info.first.second, table_info.second);
    }

    for (auto table_info : leaf_edges) {
        printf("LeafTable %u-%u: %d\n", table_info.first.first, table_info.first.second, table_info.second);
    }
#ifdef SPARSE_BITMAP
    for (auto table_info : core_edges_compact_ratio) {
        printf("CoreTable %u-%u: %.6f, %u\n", table_info.first.first, table_info.first.second, table_info.second.first, table_info.second.second);
    }
#endif
    printf("Total Cardinality: %.1lf\n", sum);
}

void catalog::initialize_catalog_info() {
     for (uint32_t u = 0; u < num_sets_; ++u) {
        for (uint32_t v = u + 1; v < num_sets_; ++v) {
            if (query_graph_->checkEdgeExistence(u, v)) {
                RelationMetaInfo meta_info;
                if (query_graph_->getCoreValue(u) > 1 && query_graph_->getCoreValue(v) > 1) {
                    meta_info.type = EdgeType::CoreEdge;
                }
                else if (query_graph_->getVertexDegree(u) == 1 || query_graph_->getVertexDegree(v) == 1) {
                    meta_info.type = EdgeType::LeafEdge;
                }
                else {
                    meta_info.type = EdgeType::TreeEdge;
                }
                catalog_info_.insert(std::make_pair(std::make_pair(u, v), meta_info));
            }
        }
    }
}
