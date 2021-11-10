//
// Created by ssunah on 6/22/18.
//

#ifndef SUBGRAPHMATCHING_CONFIG_H
#define SUBGRAPHMATCHING_CONFIG_H

/**
 * Setting the value as 1 is to (1) enable the neighbor label frequency filter (i.e., NLF filter); and (2) enable
 * to check the existence of an edge with the label information. The cost is to (1) build an unordered_map for each
 * vertex to store the frequency of the labels of its neighbor; and (2) build the label neighbor offset.
 * If the memory can hold the extra memory cost, then enable this feature to boost the performance. Otherwise, disable
 * it by setting this value as 0.
 */
#define OPTIMIZED_LABELED_GRAPH 1

#define MAXIMUM_QUERY_GRAPH_SIZE 32
#define HASH_TABLE_RATIO 1.2

/**
 * Set intersection method.
 * 0: Hybrid method; 1: Merge based set intersections.
 */
#define HYBRID 0

/**
 * Accelerate set intersection with SIMD instructions.
 * 0: AVX2; 1: AVX512; 2: Basic;
 */
#define SI 0

/**
 * Optimize relation structure.
 * 0: encoded trie; 1: hash; 2: trie.
 */
#define RELATION_STRUCTURE 0

/**
 * Add the macro to find homomorphisms.
 */
// #define HOMOMORPHISM

/**
 * Add the macro to enable intersection caching.
 */
// #define INTERSECTION_CACHE

/**
 * Add the macro to enable failing set pruning.
 */
#define FAILING_SET_PRUNING

/**
 * Add the macro to enable QFilter.
 */
// #define SPARSE_BITMAP

#define COLLECT_FAIL_SI_STATISTICS
#define COLLECT_INVALID_PR_STATISTICS
#define COLLECT_SI_STATISTICS
#define OUTPUT_OPTIMIZATION

#define NANOSECTOSEC(elapsed_time) ((elapsed_time)/(double)1000000000)

// #define ENABLE_OUTPUT
#define OUTPUT_RESULT_NUM_LIMIT (1024 * 1024)

#endif //SUBGRAPHMATCHING_CONFIG_H
