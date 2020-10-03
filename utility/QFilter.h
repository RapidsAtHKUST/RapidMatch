#ifndef SUBGRAPHMATCHING_QFLITER_H
#define SUBGRAPHMATCHING_QFLITER_H

#pragma once

#include "configuration/config.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include <omp.h> //OpenMP

#include <iostream>
#include <vector>
#include <fstream>
#include <set>
#include <atomic> //CAS
#include <chrono>

#include "han/intersection_algos.hpp"
#include "sparsepp/spp.h"

using spp::sparse_hash_map;
using namespace std;
using namespace std::chrono;

struct BSRGraph {
    vector<BSRSet> bsrs;
    sparse_hash_map<uint32_t, BSRSet> hash_bsrs_;
    int max_d_ = 0;

    BSRGraph() = default;
    ~BSRGraph() {
        for (auto iter : bsrs) {
            delete[] iter.base_;
            delete[] iter.states_;
        }

        for (auto iter : hash_bsrs_) {
            delete[] iter.second.base_;
            delete[] iter.second.states_;
        }
    }

    template<typename OFF, typename T>
    void load(size_t num_vertices, OFF &off_beg, OFF &node_off_end, T &adj, uint32_t max_degree, bool hashed) {
        bsrs.resize(num_vertices);

        max_d_ = max_degree;
        {
            if (node_off_end != off_beg) {
                cout << "err" << endl;
                for (auto u = 0u; u < num_vertices; u++) {
                    node_off_end[u + 1] = static_cast<uint32_t>(
                            lower_bound(adj + off_beg[u], adj + off_beg[u + 1], u) - adj);
                }
            }

            int *tmp_base = new int[max_d_];
            int *tmp_state = new int[max_d_];

            for (int i = 0; i < num_vertices; i++) {

                auto degree = node_off_end[i + 1] - off_beg[i];
                auto tmp_size = offline_uint_trans_bsr(reinterpret_cast<int *>(adj) + off_beg[i], degree, tmp_base,
                                                       tmp_state);

                bsrs[i].base_ = new int[degree + 16];
                bsrs[i].states_ = new int[degree + 16];
                bsrs[i].size_ = tmp_size;
                memcpy(bsrs[i].base_, tmp_base, static_cast<size_t>(tmp_size) * sizeof(int));
                memcpy(bsrs[i].states_, tmp_state, static_cast<size_t>(tmp_size) * sizeof(int));
            }

            delete[] tmp_base;
            delete[] tmp_state;
        }
    }

    void load(uint32_t num_vertices, uint32_t *parents, uint32_t *nbr_offset, uint32_t *nbrs, uint32_t max_degree) {
        hash_bsrs_.reserve(num_vertices);
        max_d_ = max_degree;
        {
            int* tmp_base = new int[max_degree];
            int* tmp_state = new int[max_degree];

            for (uint32_t i = 0; i < num_vertices; ++i) {
                uint32_t key = parents[i];
                uint32_t degree = nbr_offset[i + 1] - nbr_offset[i];
                uint32_t offset = nbr_offset[i];
                auto tmp_size = offline_uint_trans_bsr(reinterpret_cast<int *>(nbrs) + offset, degree, tmp_base,
                                                       tmp_state);

                BSRSet temp_bsr;
                temp_bsr.base_ = new int[degree + 16];
                temp_bsr.states_ = new int[degree + 16];
                temp_bsr.size_ = tmp_size;
                memcpy(temp_bsr.base_, tmp_base, static_cast<size_t>(tmp_size) * sizeof(int));
                memcpy(temp_bsr.states_, tmp_state, static_cast<size_t>(tmp_size) * sizeof(int));
                hash_bsrs_.emplace(key, temp_bsr);
            }

            delete[] tmp_base;
            delete[] tmp_state;
        }
    }

    void load(sparse_hash_map<uint32_t, std::pair<uint32_t, uint32_t>>* adj, uint32_t* nbrs, uint32_t max_degree) {
        hash_bsrs_.reserve(adj->size());

        max_d_ = max_degree;
        {
            int* tmp_base = new int[max_degree];
            int* tmp_state = new int[max_degree];

            for (auto& iter : *adj) {
                auto key = iter.first;
                auto degree = iter.second.first;
                auto offset = iter.second.second;
                auto tmp_size = offline_uint_trans_bsr(reinterpret_cast<int *>(nbrs) + offset, degree, tmp_base,
                                                       tmp_state);

                BSRSet temp_bsr;
                temp_bsr.base_ = new int[degree + 16];
                temp_bsr.states_ = new int[degree + 16];
                temp_bsr.size_ = tmp_size;
                memcpy(temp_bsr.base_, tmp_base, static_cast<size_t>(tmp_size) * sizeof(int));
                memcpy(temp_bsr.states_, tmp_state, static_cast<size_t>(tmp_size) * sizeof(int));
                hash_bsrs_.emplace(key, temp_bsr);
            }
            delete[] tmp_base;
            delete[] tmp_state;
        }
    }

};

#endif //SUBGRAPHMATCHING_QFLITER_H
