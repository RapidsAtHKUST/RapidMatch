
#include "nop_join.h"
#include <utility>



    uint64_t nop_join::execute() {
        uint64_t size_l = left_->get_size();
        uint64_t size_r = right_->get_size();

        if(size_l == 0 || size_r == 0){
            return 0;
        }

        if (size_l > size_r) {
            swapped_ = true;
            std::swap(size_l, size_r);
            std::swap(left_, right_);
            std::swap(lkp1_, rkp1_);
            std::swap(lkp2_, rkp2_);
        }

        uint32_t l_length = left_->get_tuple_length();
        uint32_t r_length = right_->get_tuple_length();

        uint64_t result_count = 0;
        auto new_size = static_cast<uint64_t>(hash_table_size_ratio_ * size_l);
        hash_table table = hash_table(new_size);
        // Build Phase
        for(uint64_t k = 0; k < size_l; ++k){
            uint32_t* data = left_->get_tuple(k);

            uint64_t key;
            uint64_t index;
            if (kc_ == 1) {
                uint32_t temp_key = get_key(data, lkp1_);
                index = murmur3(temp_key) % new_size;
                key = temp_key;
            }
            else {
                uint64_t temp_key = get_key(data, lkp1_, lkp2_);
                index = murmur3(temp_key) % new_size;
                key = temp_key;
            }

            my_tuple curr = std::make_tuple(key, k);
            hash_table::bucket& bucket = table.arr[index];
            switch(bucket.count){
                case 0:
                    bucket.t1 = curr;
                    break;
                case 1:
                    bucket.t2 = curr;
                    break;
                case 2:
                    bucket.next = std::make_unique<overflow>(curr);
                    break;
                default:
                    overflow* ptr = bucket.next.get();
                    // Follow pointer indirection
                    uint32_t overflow_tuple_count = bucket.count - 2;
                    uint32_t overflow_bucket_count = (overflow_tuple_count - 1) / TUPLE_COUNT_PER_OVERFLOW_BUCKET;
                    uint32_t overflow_bucket_offset = overflow_tuple_count % TUPLE_COUNT_PER_OVERFLOW_BUCKET;

                    for(uint64_t i = 0; i < overflow_bucket_count; i++){
                        ptr = ptr->next.get();
                    }

                    if (overflow_bucket_offset != 0) {
                        // The last bucket still has space.
                        ptr->tuples[overflow_bucket_offset] = curr;
                    }
                    else {
                        ptr->next = std::make_unique<overflow>(curr);
                    }

            }
            ++bucket.count;
        }
        // Probe Phase
        for(uint64_t k = 0; k < size_r; ++k){
            uint32_t* data = right_->get_tuple(k);

            uint64_t key;
            uint64_t index;
            if (kc_ == 1) {
                uint32_t temp_key = get_key(data, rkp1_);
                index = murmur3(temp_key) % new_size;
                key = temp_key;
            }
            else {
                uint64_t temp_key = get_key(data, rkp1_, rkp2_);
                index = murmur3(temp_key) % new_size;
                key = temp_key;
            }

            hash_table::bucket& bucket = table.arr[index];

            // Look at first tuple
            if(bucket.count > 0 && std::get<0>(bucket.t1) == key){
                result_count += 1;
                if (materialize_) {
                    uint64_t l_offset = std::get<1>(bucket.t1);
                    uint64_t r_offset = k;
                    if (kc_ == 1) {
                        if (!swapped_) {
                            reconstruct(left_->get_tuple(l_offset), l_length,
                                        right_->get_tuple(r_offset), r_length, rkp1_);
                        }
                        else {
                            reconstruct(right_->get_tuple(r_offset), r_length, left_->get_tuple(l_offset),
                                        l_length, lkp1_);
                        }
                    } else {
                        if (!swapped_) {
                            reconstruct(left_->get_tuple(l_offset), l_length,
                                        right_->get_tuple(r_offset), r_length, rkp1_, rkp2_);
                        }
                        else {
                            reconstruct(right_->get_tuple(l_offset), r_length,
                                        left_->get_tuple(r_offset), l_length, lkp1_, lkp2_);
                        }

                    }
                }
            }

            // Look at second tuple
            if(bucket.count > 1 && std::get<0>(bucket.t2) == key){
                result_count += 1;
                if (materialize_) {
                    uint64_t l_offset = std::get<1>(bucket.t2);
                    uint64_t r_offset = k;
                    if (kc_ == 1) {
                        if (!swapped_) {
                            reconstruct(left_->get_tuple(l_offset), l_length,
                                        right_->get_tuple(r_offset), r_length, rkp1_);
                        }
                        else {
                            reconstruct(right_->get_tuple(r_offset), r_length, left_->get_tuple(l_offset),
                                        l_length, lkp1_);
                        }
                    } else {
                        if (!swapped_) {
                            reconstruct(left_->get_tuple(l_offset), l_length,
                                        right_->get_tuple(r_offset), r_length, rkp1_, rkp2_);
                        }
                        else {
                            reconstruct(right_->get_tuple(l_offset), r_length,
                                        left_->get_tuple(r_offset), l_length, lkp1_, lkp2_);
                        }

                    }

                }
            }

            // Follow overflow buckets
            if(bucket.count > 2){
                overflow* curr_over = bucket.next.get();
                for(uint32_t i = 0; i < bucket.count - 2; ++i){
                    uint32_t overflow_bucket_offset = i % TUPLE_COUNT_PER_OVERFLOW_BUCKET;
                    if(std::get<0>(curr_over->tuples[overflow_bucket_offset]) == key){
                        result_count += 1;
                        if (materialize_) {
                            uint64_t l_offset = std::get<1>(curr_over->tuples[overflow_bucket_offset]);
                            uint64_t r_offset = k;
                            if (kc_ == 1) {
                                if (!swapped_) {
                                    reconstruct(left_->get_tuple(l_offset), l_length,
                                                right_->get_tuple(r_offset), r_length, rkp1_);
                                }
                                else {
                                    reconstruct(right_->get_tuple(r_offset), r_length, left_->get_tuple(l_offset),
                                                l_length, lkp1_);
                                }
                            } else {
                                if (!swapped_) {
                                    reconstruct(left_->get_tuple(l_offset), l_length,
                                                right_->get_tuple(r_offset), r_length, rkp1_, rkp2_);
                                }
                                else {
                                    reconstruct(right_->get_tuple(l_offset), r_length,
                                                left_->get_tuple(r_offset), l_length, lkp1_, lkp2_);
                                }

                            }
                        }
                    }
                    if (overflow_bucket_offset== TUPLE_COUNT_PER_OVERFLOW_BUCKET - 1) {
                        curr_over = curr_over->next.get();
                    }
                }
            }
        }
        return result_count;
    }

nop_join::nop_join(flat_relation *left, uint32_t lkp, flat_relation *right, uint32_t rkp, flat_relation *output,
                   double table_size) : left_(left), lkp1_(lkp), right_(right), rkp1_(rkp),
                                        output_(output), hash_table_size_ratio_(table_size), kc_(1), swapped_(false), materialize_(false){
    materialize_ = output != nullptr;
}

nop_join::nop_join(flat_relation *left, uint32_t lkp1, uint32_t lkp2, flat_relation *right, uint32_t rkp1, uint32_t rkp2,
                   flat_relation *output, double table_size) : left_(left), lkp1_(lkp1), lkp2_(lkp2),
                                                               right_(right), rkp1_(rkp1), rkp2_(rkp2),
                                                               output_(output), hash_table_size_ratio_(table_size), kc_(2), swapped_(false), materialize_(false){
    materialize_ = output != nullptr;
}


