/*
 * Written by Scott Norton
 * 15-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#ifndef DRUGSEQ_DEDUPLICATE_DUPLICATEDETECTOR_H
#define DRUGSEQ_DEDUPLICATE_DUPLICATEDETECTOR_H

#include <string>
#include <unordered_map>
#include <optional>
#include "FilterKey.hpp"

class DuplicateDetector {
    using memo_pos_t = std::unordered_map<std::string, std::string>;
    using memo_pos_pair_t = std::pair<const std::string, std::string>;
    using memo_chrom_t = std::unordered_map<unsigned long long, memo_pos_t>;
    using memo_chrom_pair_t = std::pair<const unsigned long long, memo_pos_t>;
    using memo_t = std::unordered_map<ChromKey, memo_chrom_t>;
    using memo_pair_t = std::pair<const ChromKey, memo_chrom_t>;

    // From constructor
    const unsigned long long& umilen;
    const unsigned long long& max_distance;
    const unsigned long long& max_lookback;

    memo_t memo; // UMI memory
    std::vector<double> nw_memory; // for NW algorithm
    double min_score; // minimum alignment score between two UMIs to consider them duplicates. Computed as umilen - max_distance.

    // Implements the Needleman-Wunsch algorithm
    double aligner_score(const std::string& a, const std::string& b);
public:
    // Class constructor
    //   const int& umilen - Fixed length of the UMI, in bases
    //   const int& max_lookback - Width of alignment error window
    //   const int& max_distance - Maximum distance between two UMIs to consider them duplicates
    DuplicateDetector(const unsigned long long &umilen, const unsigned long long &max_lookback, const unsigned long long &max_distance);

    // Clear all UMIs from memory
    void reset();

    // Check whether the read encoded in key is a new UMI
    //   const FilterKey& key - Struct with details of the read's alignments and UMI sequence
    // returns: std::optional<std::string> - If key is a duplicate read, the query name of the corresponding original read.
    std::optional<std::string> check_umi(const FilterKey& key);
};

#endif //DRUGSEQ_DEDUPLICATE_DUPLICATEDETECTOR_H
