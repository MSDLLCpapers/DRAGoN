/*
 * Written by Scott Norton
 * 12-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <bamtools/api/BamConstants.h>
#include "BarcodeMatcher.hpp"
#include "SeqToMask.hpp"

// Compare two binarized DNA sequences with at most this->max_mismatch substitutions.
//  Args:
//    mask: Left-hand side
//    refmask: Right-hand side
//  Returns:
//    Number of substitutions. Note that this function will return the moment it detects that the number of substitutions will exceed this->max_mismatch, so this will be a lower bound on the actual number of mismatches between the two sequences.
//  Comment: The parameter names refer to the fact that the binarization effectively turns each base into a bitmask where each bit represents one of the canonical DNA bases (A, C, G, T). The ambiguity code R has bits 0 and 2 set, corresponding to ambiguity between A and G.
//           Thus, a bitwise AND with either base's binarized form will be nonzero.
int BarcodeMatcher::match_one(std::vector<unsigned char> const& mask, std::vector<unsigned char> const& refmask) const {
    int nmismatch = 0;
    for (int i = 0; i < bclen; ++i) {
        if (!(mask.at(i) & refmask.at(i)) && ++nmismatch > max_mismatch) {
            break;
        }
    }
    return nmismatch;
}

// Constructor for the BarcodeMatcher
//  Args:
//    filename: String containing the path to a tab-delimited file containing at least two columns. The first column is the sequence of each well's barcode, and the second is the name of that well.
//    max_mismatch: A query barcode will match a reference barcode with at most this many substitutions.
//    max_n: A query barcode with only be tested for match to a reference barcode if it has at most this many N bases.
BarcodeMatcher::BarcodeMatcher(const std::string& filename, const int &max_mismatch, const int &max_n) : max_mismatch(max_mismatch), max_n(max_n), numbc(0) {
    // Parse the file
    std::ifstream bcfile(filename);
    std::string line;

    bclen = std::string::npos;
    while (std::getline(bcfile, line)) {
        std::size_t tab1 = line.find('\t');
        if (tab1 == std::string::npos) {
            throw std::runtime_error("ERR: " + filename + ":" + std::to_string(numbc + 1) + ": not tab-delimited");
        }
        std::size_t tab2 = line.find('\t', tab1 + 1);
        if (bclen == std::string::npos) {
            bclen = tab1;
        } else if (tab1 != bclen) {
            throw std::runtime_error("ERR: " + filename + ":" + std::to_string(numbc + 1) + ": barcode length disagree");
        }
        if (tab2 == tab1 + 1) {
            throw std::runtime_error("ERR: " + filename + ":" + std::to_string(numbc + 1) + ": barcode name missing");
        }
        // Binarization
        seq_to_mask(bcseqs.emplace_back(line.substr(0, bclen)), bcmasks.emplace_back(bclen));
        bcnames.push_back(line.substr(bclen + 1, tab2 == std::string::npos ? tab2 : tab2 - bclen - 1));
        ++numbc;
    }
    nobcseq.assign(bclen, 'N');
    // ~bcfile
}

// Search for the first reference barcode closely matching the query sequence
//  Args:
//    seq: The query DNA sequence
//  Returns:
//    A pair<int, int> containing the index of the matched barcode and the number of substitutions between it and the query.
std::pair<int, int> const& BarcodeMatcher::match_barcode(std::string const& seq) {
    auto it = memo.find(seq);
    if (it == memo.cend()) {
        std::pair<int, int> &result = memo[seq];
        result = std::make_pair(-1, -1);

        std::vector<unsigned char> mask(bclen);
        if (std::count(seq.cbegin(), seq.cend(), 'N') <= max_n) {
            // Test the query barcode against each sequence in the reference
            seq_to_mask(seq, mask);
            int best_score = max_mismatch + 1;
            for (auto i_refmask = bcmasks.begin(); i_refmask != bcmasks.end(); ++i_refmask) {
                int score = match_one(mask, *i_refmask);
                if (score < best_score) {
                    result = std::make_pair(static_cast<int>(i_refmask - bcmasks.begin()), best_score = score);
                }
            }
        }
    }
    return memo.at(seq);
}
