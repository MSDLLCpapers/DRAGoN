/*
 * Written by Scott Norton
 * 12-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#ifndef DRUGSEQ_QUALITYFILTER_H
#define DRUGSEQ_QUALITYFILTER_H

#include <string>
#include <array>
#include <unordered_set>
#include <optional>
#include <vector>
#include "ReadPair.hpp"
#include "QCStats.hpp"

struct CutAdaptDPMatrixCell {
    int cost;
    int score;    // score for this alignment (mostly keeps track of matches)
    int origin;   // where the alignment originated: negative for positions within seq1, positive for pos. within seq2

    // for use with std::iota
    CutAdaptDPMatrixCell& operator++() {
        ++cost;
        return *this;
    }
};

class QCFiltrator {
    bool enable_length = false;
    bool enable_umiqual = false;
    bool enable_homogeny = false;
    bool enable_homopolymer = false;
    bool enable_adapter_trimming = false;
    bool enable_trimming = false;
    bool enable_duplicate_filter = false;

    int min_len;
    int min_qual;
    int umi_qual_bases;
    int homop_len;
    std::string adapter;
    double max_AT, max_GC;
    std::array<std::string, 5> homopolymers;
    std::unordered_set<ReadPair> seen;
    std::vector<CutAdaptDPMatrixCell> adapter_alignment_matrix;
    std::optional<CutAdaptDPMatrixCell> cutadapt_match_adapter(ReadPair const& rp);
public:
    QCFiltrator(const int &min_len, const int &min_qual, const int &phred, const int &umi_qual_bases, const int &homop_len, const double &max_AT, const double &max_GC, const std::string &adapter);
    QCFailReason operator()(ReadPair& rp);
};

#endif //DRUGSEQ_QUALITYFILTER_H
