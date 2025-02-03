/*
 * Written by Scott Norton
 * 15-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#include <string>
#include <unordered_map>
#include <algorithm>  // std::max(initializer-list)
#include <optional>
#include "FilterKey.hpp"
#include "DuplicateDetector.hpp"

DuplicateDetector::DuplicateDetector(const unsigned long long &umilen, const unsigned long long &max_distance, const unsigned long long &max_lookback) :
        umilen(umilen),
        max_lookback(max_lookback),
        max_distance(max_distance),
        min_score(umilen - max_distance),
        nw_memory(umilen * 2, 0)
{}

void DuplicateDetector::reset() {
    memo.clear();
}

double DuplicateDetector::aligner_score(const std::string& a, const std::string& b) {
    // Needleman-Wunsch algorithm
    // In this implementation, all gap penalties are the same, whether open, extend, or close, left or right.

    constexpr double s_match = 1.0;
    constexpr double s_mismatch = 0.0;
    constexpr double s_gap = -0.5;

    // To keep memory requirements low, we use two halves of a vector of size umilen
    // and swap which one is the current row and which is the previous in the
    // dynamic programming array.
    auto iter_a = nw_memory.begin(), iter_b = iter_a + umilen;
    auto p_iter_prev = &iter_a, p_iter_curr = &iter_b;

    // Partially unroll the loop
    nw_memory[0] = a[0] == b[0] ? s_match : s_mismatch;
    for (size_t j = 1; j < umilen; ++j) {
        nw_memory[j] = std::max({
            nw_memory[j - 1] + s_gap,
            s_gap * j + (a[0] == b[j] ? s_match : s_mismatch)
        });
    }

    for (size_t i = 1; i < umilen; ++i) {
        auto iter_prev = *p_iter_prev, iter_curr = *p_iter_curr, iter_prev_m1 = iter_prev, iter_curr_m1 = iter_curr;
        *iter_curr++ = std::max({
            *iter_prev++ + s_gap,
            s_gap * i + (a[i] == b[0] ? s_match : s_mismatch)
        });
        for (size_t j = 1; j < umilen; ++j) {
            *iter_curr = std::max({
                *iter_prev + s_gap,
                *iter_curr_m1 + s_gap,
                *iter_prev_m1 + (a[i] == b[j] ? s_match : s_mismatch)
            });
            iter_curr_m1 = iter_curr++;
            iter_prev_m1 = iter_prev++;
        }
        std::swap(p_iter_curr, p_iter_prev);
    }

    return *(*p_iter_prev + umilen - 1);
}

std::optional<std::string> DuplicateDetector::check_umi(FilterKey const& key) {
    for (const FilterKey::coord_pair_t &chrom : key.coords) {
        // If we've not yet seen a read from this contig/strand, skip and move on.
        memo_t::const_iterator chrom_memo_it = memo.find(chrom.first);
        if (chrom_memo_it == memo.cend()) {
            continue;
        }
        // Survey every start position this read maps to on this contig/strand.
        for (const unsigned long long &pos : chrom.second) {
            // Look +/-(max_lookback) from the read position.
            for (unsigned long long i = pos > max_lookback ? pos - max_lookback : 0; i <= pos + max_lookback; ++i) {
                // If we've not yet seen a read from this coordinate, skip and move on.
                memo_chrom_t::const_iterator pos_memo_it = chrom_memo_it->second.find(i);
                if (pos_memo_it == chrom_memo_it->second.cend()) {
                    continue;
                }
                // Check for exact match
                memo_pos_t::const_iterator umi_memo_it = pos_memo_it->second.find(key.umi);
                if (umi_memo_it != pos_memo_it->second.cend()) {
                    return umi_memo_it->second;
                }
                // If max_distance is 0, restrict to exact matches only.
                if (max_distance == 0) {
                    continue;
                }
                // Check for close matches
                for (const std::pair<const std::string, std::string> &umi_pair : pos_memo_it->second) {
                    if (aligner_score(key.umi, umi_pair.first) >= min_score) {
                        return umi_pair.second;
                    }
                }
            }
        }
    }
    // New UMI, record it in memory and return the unique sentinel
    for (const FilterKey::coord_pair_t &chrom : key.coords) {
        memo_chrom_t &chrom_memo = memo[chrom.first];
        for (const unsigned long long &pos : chrom.second) {
            chrom_memo[pos][key.umi] = key.q_name;
        }
    }
    return {};
}
