#ifndef SUBGRAPHMATCHING_SEARCH_H
#define SUBGRAPHMATCHING_SEARCH_H

#include <cstdint>

class search {
public:
    static uint32_t binary_search(const uint32_t *src, uint32_t begin, uint32_t end, uint32_t target);
};


#endif //SUBGRAPHMATCHING_SEARCH_H
