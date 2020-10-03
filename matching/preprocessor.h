#ifndef SUBGRAPHMATCHING_PREPROCESSOR_H
#define SUBGRAPHMATCHING_PREPROCESSOR_H


#include <cstdint>
#include "graph/graph.h"
#include "relation/catalog.h"

class preprocessor {
public:
    double preprocess_time_;
    double filter_time_;
    double scan_time_;
    double semi_join_time_;
    double bottom_up_non_core_semi_join_time_;
    double bottom_up_core_semi_join_time_;
    double top_down_non_core_semi_join_time_;
    double top_down_core_semi_join_time_;

private:
    uint32_t vertices_count_;
    uint32_t non_core_vertices_count_;
    uint32_t* degeneracy_ordering_;
    uint32_t* vertices_index_;
    uint32_t* non_core_vertices_parent_;
    uint32_t* non_core_vertices_children_;
    uint32_t* non_core_vertices_children_offset_;

    const Graph* query_graph_;
    const Graph* data_graph_;

private:
    void initialize(const Graph* query_graph, const Graph* data_graph);

    void scan_relation(catalog* storage);

    void eliminate_dangling_tuples(catalog *storage);

    void generate_preprocess_plan();

    edge_relation* get_key_position_in_relation(uint32_t u, uint32_t v, catalog *storage, uint32_t &kp);

public:
    preprocessor() : vertices_count_(0), non_core_vertices_count_(0), degeneracy_ordering_(nullptr), vertices_index_(nullptr),
                     non_core_vertices_parent_(nullptr), non_core_vertices_children_(nullptr), non_core_vertices_children_offset_(nullptr) {}
    ~preprocessor() {
        delete[] degeneracy_ordering_;
        delete[] vertices_index_;
        delete[] non_core_vertices_parent_;
        delete[] non_core_vertices_children_;
        delete[] non_core_vertices_children_offset_;
    }

    void execute(const Graph *query_graph, const Graph *data_graph, catalog *storage, bool enable_elimination);

    void print_metrics();
};


#endif //SUBGRAPHMATCHING_PREPROCESSOR_H
