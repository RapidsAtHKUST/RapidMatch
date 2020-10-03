#ifndef SUBGRAPHMATCHING_PROJECTION_H
#define SUBGRAPHMATCHING_PROJECTION_H


#include <cstdint>
#include <relation/edge_relation.h>
/// Will sort the result after projection.
class projection {
private:
    bool* index_;
    uint32_t* buffer_;
    uint32_t index_size_;

public:
    projection(uint32_t index_size) {
        index_size_ = index_size;
        buffer_ = new uint32_t[index_size_];
        index_ = new bool[index_size_];
    }
    ~projection() {
        delete[] buffer_;
        delete[] index_;
    }

    void execute(edge_relation *relation, uint32_t kp, uint32_t* &res, uint32_t& res_cnt);
};


#endif //SUBGRAPHMATCHING_PROJECTION_H
