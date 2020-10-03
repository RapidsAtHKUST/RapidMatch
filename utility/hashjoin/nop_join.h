//
// Benjamin Wagner 2018
//

#ifndef HASHJOINS_NOP_JOIN_H
#define HASHJOINS_NOP_JOIN_H

#include <memory>
#include <vector>
#include <tuple>
#include "../relation/flat_relation.h"
#include "hash_helpers.h"
#include "../pretty_print.h"
    /// Single threaded no partitioning join
class nop_join {
private:
    flat_relation* left_;
    uint32_t lkp1_;
    uint32_t lkp2_;
    flat_relation* right_;
    uint32_t rkp1_;
    uint32_t rkp2_;
    flat_relation* output_;
    double hash_table_size_ratio_;
    uint32_t kc_;
    bool swapped_;
    bool materialize_;

public:
    nop_join(flat_relation *left, uint32_t lkp, flat_relation *right, uint32_t rkp, flat_relation *output,
                         double table_size);

    nop_join(flat_relation *left, uint32_t lkp1, uint32_t lkp2, flat_relation *right, uint32_t rkp1, uint32_t rkp2,
                 flat_relation *output, double table_size);

    uint64_t execute();

private:
    inline uint64_t get_key(uint32_t* data, uint32_t kp1, uint32_t kp2) {
        uint64_t key = (uint64_t) data[kp1] << 32 | data[kp2];
        return key;
    }

    inline uint32_t get_key(uint32_t* data, uint32_t kp1) {
        return data[kp1];
    }

    inline void reconstruct(uint32_t* ldata, uint32_t llength, uint32_t* rdata, uint32_t rlength, uint32_t rkp) {
        auto buffer = output_->push_back();

        uint32_t offset = 0;
        memcpy(buffer + offset, ldata, sizeof(uint32_t) * llength);
        offset += llength;
        memcpy(buffer + offset, rdata, sizeof(uint32_t) * rkp);
        offset += rkp;
        memcpy(buffer + offset, rdata + rkp + 1, sizeof(uint32_t) * (rlength - rkp - 1));

//        std::cout << "Left: " << pretty_print_array(ldata, llength)
//                  << ", Right: " << pretty_print_array(rdata, rlength) << ", Result:" << pretty_print_array(buffer, llength + rlength - 1) << std::endl;
    }

    inline void reconstruct(uint32_t* ldata, uint32_t llength, uint32_t* rdata, uint32_t rlength, uint32_t rkp1, uint32_t rkp2) {
            auto buffer = output_->push_back();

            uint32_t offset = 0;
            memcpy(buffer + offset, ldata, sizeof(uint32_t) * llength);
            offset += llength;
            memcpy(buffer + offset, rdata, sizeof(uint32_t) * rkp1);
            offset += rkp1;
            memcpy(buffer + offset, rdata + rkp1 + 1, sizeof(uint32_t) * (rkp2 - rkp1 - 1));
            offset += rkp2 - rkp1 - 1;
            memcpy(buffer + offset, rdata + rkp2 + 1, sizeof(uint32_t) * (rlength - rkp2 - 1));

//            std::cout << "Left: " << pretty_print_array(ldata, llength)
//                      << ", Right: " << pretty_print_array(rdata, rlength) << ", Result:" << pretty_print_array(buffer, llength + rlength - 2) << std::endl;
    }

};

#endif  // HASHJOINS_NOP_JOIN_H
