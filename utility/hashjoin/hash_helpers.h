//
// Benjamin Wagner 2018
//

#ifndef HASHJOINS_HASH_HELPERS_H
#define HASHJOINS_HASH_HELPERS_H

#include <memory>
#include <tuple>
#include <mutex>
#include <vector>

#define TUPLE_COUNT_PER_OVERFLOW_BUCKET 4
/// Murmur 3 hash function for 64 bit uint keys, inline to keep all in header file
inline uint64_t murmur3(uint64_t val) {
        // Murmur 3 taken from "A Seven-Dimensional Analysis of Hashing Methods and its
        // Implications on Query Processing" by Richter et al
        val ^= val >> 33;
        val *= 0xff51afd7ed558ccd;
        val ^= val >> 33;
        val *= 0xc4ceb9fe1a85ec53;
        val ^= val >> 33;
        return val;
}

inline uint32_t murmur3(uint32_t val) {
    val ^= val >> 16;
    val *= 0x85ebca6b;
    val ^= val >> 13;
    val *= 0xc2b2ae35;
    val ^= val >> 16;
    return val;
}


/// First element is the value on which should be joined, second one the rid
typedef std::tuple<uint64_t, uint64_t> my_tuple;

/// Overflow bucket used for chaining, we set the size of overflow bucket as 16.
struct overflow {
    my_tuple tuples[TUPLE_COUNT_PER_OVERFLOW_BUCKET];
    std::unique_ptr<overflow> next;
    explicit overflow(my_tuple t) { tuples[0] = t; }
};

/// Simple chained hash table used within the nop join
struct hash_table{
    /// One of the hash table entries
    struct bucket{
        uint32_t count;
        my_tuple t1;
        my_tuple t2;
        std::unique_ptr<overflow> next;

        /// Default constructor
        bucket(): count(0), next(nullptr) {}
    };

    std::unique_ptr<bucket[]> arr;
    uint64_t size;

    explicit hash_table(uint64_t size): size(size){
        arr = std::make_unique<bucket[]>(size);
    }
};

#endif  // HASHJOINS_HASH_HELPERS_H
