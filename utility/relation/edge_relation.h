#ifndef SUBGRAPHMATCHING_EDGE_RELATION_H
#define SUBGRAPHMATCHING_EDGE_RELATION_H

#include <cstdint>
#include "configuration/types.h"


class edge_relation {
public:
    edge* edges_;
    uint32_t size_;

    edge_relation() : edges_(nullptr), size_(0) {}
    ~edge_relation() {
        delete [] edges_;
    }

    uint64_t memory_cost() {
        return sizeof(edge) * size_;
    }
};

#endif //SUBGRAPHMATCHING_EDGE_RELATION_H
