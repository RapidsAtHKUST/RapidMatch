//
// Created by ssunah on 12/19/19.
//

#include <vector>
#include "nop_join.h"
#include "../pretty_print.h"
int main() {
    {
        // Test 1 join key.
        const uint32_t unique_key_count = 3;
        const uint32_t per_unique_key_count = 16;

        flat_relation left(3, 65536);
        flat_relation right(3, 65536);
        flat_relation output(5, 65536);

        std::vector<std::vector<uint32_t>> left_answer;
        std::vector<std::vector<uint32_t>> right_answer;
        uint32_t lkp = 0;
        uint32_t rkp = 1;
        // Generate tuples in left relation.
        for (uint32_t i = 0; i < unique_key_count; ++i) {
            for (uint32_t j = 0; j < per_unique_key_count; ++j) {
                std::vector<uint32_t> tuple(3);
                tuple[0] = i;
                tuple[1] = j;
                tuple[2] = i + j;
                left.push_back(tuple.data());
                left_answer.emplace_back(tuple);
            }
        }

        // Generate tuples in right relation.
        for (uint32_t i = 0; i < unique_key_count; ++i) {
            for (uint32_t j = 0; j < per_unique_key_count; ++j) {
                std::vector<uint32_t> tuple(3);
                tuple[0] = j + 10;
                tuple[1] = i;
                tuple[2] = i + j + 10;
                right.push_back(tuple.data());
                right_answer.emplace_back(tuple);
            }
        }

        // Join.
        std::vector<std::vector<uint32_t>> answer;
        for (auto& r_tuple : right_answer) {
            for (auto& l_tuple : left_answer) {
                if (l_tuple[lkp] == r_tuple[rkp]) {
                    answer.emplace_back(l_tuple);
                    answer.back().push_back(r_tuple[0]);
                    answer.back().push_back(r_tuple[2]);
                }
            }
        }

        nop_join join(&left, lkp, &right, rkp, &output, 1.2);
        uint64_t res_count = join.execute();

        std::cout << "Expected Number of Results: " << answer.size()
                  << " - Actual Number of Results: " << output.get_tuple_count() << std::endl;

        assert(res_count == output.get_tuple_count());
        assert(answer.size() == output.get_tuple_count());

        for (uint64_t i = 0; i < res_count; ++i) {
            auto& expected = answer[i];
            auto* value = (uint32_t*)output.get_tuple(i);

            for (uint64_t j = 0; j < 5; ++j) {
                // std::cout << "(" << i << ',' << j << ")-" << expected[j] << '-' << value[j] << std::endl;
                assert(expected[j] == value[j]);
            }
        }
    }

    {
        // Test 2 join keys.
        const uint32_t unique_key_count = 3;
        const uint32_t per_unique_key_count = 16;

        flat_relation left(5, 65536);
        flat_relation right(5, 65536);
        flat_relation output(8, 65536);

        std::vector<std::vector<uint32_t>> left_answer;
        std::vector<std::vector<uint32_t>> right_answer;
        uint32_t lkp1 = 0;
        uint32_t lkp2 = 3;
        uint32_t rkp1 = 1;
        uint32_t rkp2 = 3;
        // Generate tuples in left relation.
        for (uint32_t i = 0; i < unique_key_count; ++i) {
            for (uint32_t j = 0; j < per_unique_key_count; ++j) {
                std::vector<uint32_t> tuple(5);
                tuple[0] = i;
                tuple[1] = j;
                tuple[2] = i + j + 30;
                tuple[3] = unique_key_count - i;
                tuple[4] = i * j;
                left.push_back(tuple.data());
                left_answer.emplace_back(tuple);
            }
        }

        // Generate tuples in right relation.
        for (uint32_t i = 0; i < unique_key_count; ++i) {
            for (uint32_t j = 0; j < per_unique_key_count; ++j) {
                std::vector<uint32_t> tuple(5);
                tuple[0] = j + 10;
                tuple[1] = i;
                tuple[2] = i + j + 10;
                tuple[3] = unique_key_count - i;
                tuple[4] = i + j + 100;
                right.push_back(tuple.data());
                right_answer.emplace_back(tuple);
            }
        }

        // Join.
        std::vector<std::vector<uint32_t>> answer;
        for (auto& r_tuple : right_answer) {
            for (auto& l_tuple : left_answer) {
                if (l_tuple[lkp1] == r_tuple[rkp1] && l_tuple[lkp2] == r_tuple[rkp2]) {
                    answer.emplace_back(l_tuple);
                    answer.back().push_back(r_tuple[0]);
                    answer.back().push_back(r_tuple[2]);
                    answer.back().push_back(r_tuple[4]);
                }
            }
        }

        nop_join join(&left, lkp1, lkp2, &right, rkp1, rkp2, &output, 1.2);
        uint64_t res_count = join.execute();

        std::cout << "Expected Number of Results: " << answer.size()
                  << " - Actual Number of Results: " << output.get_tuple_count() << std::endl;

        assert(res_count == output.get_tuple_count());
        assert(answer.size() == output.get_tuple_count());

        for (uint64_t i = 0; i < res_count; ++i) {
            auto& expected = answer[i];
            auto* value = (uint32_t*)output.get_tuple(i);

            for (uint64_t j = 0; j < 8; ++j) {
                // std::cout << "(" << i << ',' << j << ")-" << expected[j] << '-' << value[j] << std::endl;
                assert(expected[j] == value[j]);
            }
        }
    }
}

