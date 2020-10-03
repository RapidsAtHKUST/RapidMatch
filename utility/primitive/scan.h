#ifndef SUBGRAPHMATCHING_SCAN_H
#define SUBGRAPHMATCHING_SCAN_H

#include "relation/edge_relation.h"
#include "graph/graph.h"
#include <vector>

class scan {
private:
    std::vector<bool> flag_;
    edge* buffer_;
    const Graph* data_graph_;

private:
    void execute_without_index(uint32_t src_label, uint32_t dst_label, edge_relation *relation);
    void execute_with_index(uint32_t src_label, uint32_t dst_label, edge_relation *relation);

public:
    scan(const Graph* data_graph) {
        data_graph_ = data_graph;
        buffer_ = new edge[data_graph_->getEdgesCount() * 2];
        flag_.resize(data_graph->getVerticesCount(), false);
    }

    ~scan() {
        delete[] buffer_;
    }

    void execute(uint32_t src_label, uint32_t dst_label, edge_relation *relation, bool indexed);

    void execute(std::vector<uint32_t>& src_candidate_set, std::vector<uint32_t>& dst_candidate_set, edge_relation *relation);
};


#endif //SUBGRAPHMATCHING_SCAN_H
