#ifndef SUBGRAPHMATCHING_SEMI_JOIN_H
#define SUBGRAPHMATCHING_SEMI_JOIN_H

#include "relation/edge_relation.h"
/// The left relation performs a semi join with the right relation.
class semi_join {
private:
    bool* index_;
    uint32_t index_size_;
public:
    semi_join(uint32_t index_size) {
        index_size_ = index_size;
        index_ = new bool[index_size];
    }
    ~semi_join() {
        delete[] index_;
    }
    uint32_t execute(edge_relation *left, const uint32_t lkp, edge_relation *right, const uint32_t rkp);
};


#endif //SUBGRAPHMATCHING_SEMI_JOIN_H
