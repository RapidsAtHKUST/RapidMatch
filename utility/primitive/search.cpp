#include "search.h"
#include <cassert>
#include <immintrin.h>

uint32_t search::binary_search(const uint32_t *src, uint32_t begin, uint32_t end, uint32_t target) {
    assert(end >= begin);
    uint32_t offset_begin = begin;
    uint32_t offset_end = end;

    while (offset_end - offset_begin >= 16) {
        auto mid = static_cast<uint32_t>((static_cast<uint64_t>(offset_begin) + offset_end) / 2);
        _mm_prefetch((char *) &src[(mid + 1 + offset_end) / 2], _MM_HINT_T0);
        _mm_prefetch((char *) &src[(mid - 1 + offset_begin) / 2], _MM_HINT_T0);
        if (src[mid] == target) {
            return mid;
        } else if (src[mid] < target) {
            offset_begin = mid + 1;
        } else {
            offset_end = mid;
        }
    }

    // linear search fallback
    for (auto offset = offset_begin; offset < offset_end; ++offset) {
        if (src[offset] > target) {
            return end;
        }
        else if (src[offset] == target) {
            return offset;
        }
    }

    return end;
}