#ifndef SUBGRAPHMATCHING_CATALOG_H
#define SUBGRAPHMATCHING_CATALOG_H

#include "encoded_trie_relation.h"
#include "hash_relation.h"
#include "trie_relation.h"
#include "edge_relation.h"
#include "../primitive/search.h"
#include "configuration/config.h"
#include "graph/graph.h"
#include "QFilter.h"
#include <map>

enum EdgeType {
    CoreEdge,
    TreeEdge,
    LeafEdge
};
struct RelationMetaInfo {
    EdgeType type;
    uint32_t after_scan_;
    uint32_t after_eliminate_;
    uint32_t max_degree_;
    uint32_t reverse_max_degree_;
};

class catalog {
public:
    std::map<std::pair<uint32_t, uint32_t>, RelationMetaInfo> catalog_info_;

public:
    uint32_t** candidate_sets_;
    uint32_t* num_candidates_;
    uint32_t num_sets_;
    uint32_t max_num_candidates_per_vertex_;
    uint32_t max_data_vertex_id_;
    Graph* query_graph_;
    Graph* data_graph_;

    hash_relation** hash_relations_;
    trie_relation** trie_relations_;
    encoded_trie_relation** encoded_trie_relations_;
    BSRGraph** bsr_relations_;
    edge_relation** edge_relations_;

private:
    void initialize_catalog_info();

public:
    explicit catalog(Graph* query_graph, Graph* data_graph) {
        query_graph_ = query_graph;
        data_graph_ = data_graph;
        num_sets_ = query_graph_->getVerticesCount();
        max_num_candidates_per_vertex_ = data_graph->getGraphMaxLabelFrequency();
        max_data_vertex_id_ = data_graph->getVerticesCount();

        candidate_sets_ = new uint32_t*[num_sets_];
        num_candidates_ = new uint32_t[num_sets_];
        edge_relations_ = new edge_relation*[num_sets_];
        hash_relations_ = new hash_relation*[num_sets_];
        trie_relations_ = new trie_relation*[num_sets_];
        encoded_trie_relations_ = new encoded_trie_relation*[num_sets_];
        bsr_relations_ = new BSRGraph*[num_sets_];

        for (uint32_t i = 0; i < num_sets_; ++i) {
            candidate_sets_[i] = nullptr;
            edge_relations_[i] = new edge_relation[num_sets_];
            hash_relations_[i] = new hash_relation[num_sets_];
            trie_relations_[i] = new trie_relation[num_sets_];
            encoded_trie_relations_[i] = new encoded_trie_relation[num_sets_];
            bsr_relations_[i] = new BSRGraph[num_sets_];
        }

        initialize_catalog_info();
    }

    ~catalog() {
        for (uint32_t i = 0; i < num_sets_; ++i) {
            delete[] candidate_sets_[i];
            delete[] edge_relations_[i];
            delete[] hash_relations_[i];
            delete[] trie_relations_[i];
            delete[] encoded_trie_relations_[i];
            delete[] bsr_relations_[i];
        }
        delete[] num_candidates_;
        delete[] candidate_sets_;
        delete[] edge_relations_;
        delete[] hash_relations_;
        delete[] trie_relations_;
        delete[] encoded_trie_relations_;
        delete[] bsr_relations_;
    }

    uint32_t get_edge_relation_cardinality(uint32_t u, uint32_t v) {
        uint32_t src = std::min(u, v);
        uint32_t dst = std::max(u, v);
        return edge_relations_[src][dst].size_;
    }

    uint32_t get_num_candidates(uint32_t u) {
        return num_candidates_[u];
    }

    uint32_t* get_candidates(uint32_t u) {
        return candidate_sets_[u];
    }

    uint32_t get_candidate(uint32_t u, uint32_t index) {
        return candidate_sets_[u][index];
    }

    uint32_t get_cardinality(uint32_t bn, uint32_t u) {
#if RELATION_STRUCTURE == 0
        return encoded_trie_relations_[bn][u].get_cardinality();
#elif RELATION_STRUCTURE == 1
        return hash_relations_[bn][u].get_cardinality();
#elif RELATION_STRUCTURE == 2
        return trie_relations_[bn][u].get_cardinality();
#endif
    }

    uint32_t get_max_degree(uint32_t bn, uint32_t u) {
#if RELATION_STRUCTURE == 0
        if (query_graph_->getCoreValue(bn) > 1 && query_graph_->getCoreValue(u) > 1)
            return encoded_trie_relations_[bn][u].max_degree();
        else
            return hash_relations_[bn][u].max_degree();
#elif RELATION_STRUCTURE == 1
        return hash_relations_[bn][u].max_degree();
#elif RELATION_STRUCTURE == 2
        return trie_relations_[bn][u].max_degree();
#endif
    }

    uint32_t get_size(uint32_t bn, uint32_t u) {
#if RELATION_STRUCTURE == 0
        return encoded_trie_relations_[bn][u].get_size();
#elif RELATION_STRUCTURE == 1
        return hash_relations_[bn][u].get_size();
#elif RELATION_STRUCTURE == 2
        return trie_relations_[bn][u].get_size();
#endif
    }

    uint32_t get_candidate_index(uint32_t u, uint32_t v) {
        return search::binary_search(candidate_sets_[u], 0, num_candidates_[u], v);
    }

    uint32_t* get_non_core_relation_children(uint32_t bn, uint32_t u, uint32_t key, uint32_t &count) {
#if RELATION_STRUCTURE == 0
        return hash_relations_[bn][u].get_children(key, count);
#elif RELATION_STRUCTURE == 1
        return hash_relations_[bn][u].get_children(key, count);
#elif RELATION_STRUCTURE == 2
        return trie_relations_[bn][u].get_children(key, count);
#endif
    }

    uint32_t* get_core_relation_children(uint32_t bn, uint32_t u, uint32_t key, uint32_t &count) {
#if RELATION_STRUCTURE == 0
        return encoded_trie_relations_[bn][u].get_children(key, count);
#elif RELATION_STRUCTURE == 1
        return hash_relations_[bn][u].get_children(key, count);
#elif RELATION_STRUCTURE == 2
        return trie_relations_[bn][u].get_children(key, count);
#endif
    }

    BSRSet get_core_relation_bsr_set(uint32_t bn, uint32_t u, uint32_t key) {
#if RELATION_STRUCTURE == 0
        return bsr_relations_[bn][u].bsrs[key];
#elif RELATION_STRUCTURE == 1
        if (!hash_relations_[bn][u].contains(key))
            return BSRSet();
        return bsr_relations_[bn][u].hash_bsrs_[key];
#elif RELATION_STRUCTURE == 2
        if (!trie_relations_[bn][u].contains(key))
            return BSRSet();
        return bsr_relations_[bn][u].hash_bsrs_[key];
#endif
    }

    bool is_adjacent(uint32_t u, uint32_t v) {
        return query_graph_->checkEdgeExistence(u, v);
    }

    uint64_t memory_cost();

    void print_cardinality();
    void print_metrics();
};

#endif //SUBGRAPHMATCHING_CATALOG_H
