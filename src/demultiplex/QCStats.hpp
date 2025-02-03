/*
 * Written by Scott Norton
 * 12-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#ifndef DRUGSEQ_QCSTATS_H
#define DRUGSEQ_QCSTATS_H

#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <nlohmann/json.hpp>
#include "ReadPair.hpp"
#include "BarcodeMatcher.hpp"

enum QCFailReason {
    PASS = 0,
    FAIL_UMI_QUALITY,
    FAIL_AT_CONTENT,
    FAIL_GC_CONTENT,
    FAIL_LENGTH,
    FAIL_HOMOPOLYMER,
    FAIL_ADAPTER_TRIMMING,
    FAIL_QUALITY_TRIMMING,
    FAIL_DUPLICATE,
    MAX,
};

static const std::string QCFailNames[QCFailReason::MAX] {
    "PASS",
    "FAIL_UMI_QUALITY",
    "FAIL_AT_CONTENT",
    "FAIL_GC_CONTENT",
    "FAIL_LENGTH",
    "FAIL_HOMOPOLYMER",
    "FAIL_ADAPTER_TRIMMING",
    "FAIL_QUALITY_TRIMMING",
    "FAIL_DUPLICATE",
};

struct Statistics {
    unsigned long long total;
    std::vector<unsigned long long> written_counts;
    std::vector<unsigned long long> well_counts;
    std::vector<unsigned long long> mismatch_counts;
    std::vector<unsigned long long> fail_counts;
    std::unordered_map<std::string, unsigned long long> unmatched_counts;
    BarcodeMatcher const& matcher;
    int max_file;
    int max_mismatch;

    Statistics(BarcodeMatcher const& matcher, const int max_file, const int &max_mismatch);
    void record(const ReadPair &rp, const QCFailReason &fail, const int &idx, const int &mismatch, const int &fileidx);
    nlohmann::ordered_json to_json() const;
    friend std::ostream &operator<<(std::ostream &strm, Statistics const &self);
};

#endif //DRUGSEQ_QCSTATS_H
