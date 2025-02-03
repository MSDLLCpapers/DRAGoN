/*
 * Written by Scott Norton
 * 12-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#ifndef DRUGSEQ_BARCODEMATCHER_H
#define DRUGSEQ_BARCODEMATCHER_H

#include <string>
#include <vector>
#include <unordered_map>

class BarcodeMatcher {
    const int &max_mismatch;
    const int &max_n;
    std::vector<std::string> bcnames;
    std::vector<std::string> bcseqs;
    std::vector<std::vector<unsigned char>> bcmasks;
    std::size_t bclen, numbc;
    std::string nobcseq;
    std::unordered_map<std::string, std::pair<int, int>> memo;

    // Compare two binarized DNA sequences with at most this->max_mismatch substitutions.
    //  Args:
    //    mask: Left-hand side
    //    refmask: Right-hand side
    //  Returns:
    //    Number of substitutions. Note that this function will return the moment it detects that the number of substitutions will exceed this->max_mismatch, so this will be a lower bound on the actual number of mismatches between the two sequences.
    //  Comment: The parameter names refer to the fact that the binarization effectively turns each base into a bitmask where each bit represents one of the canonical DNA bases (A, C, G, T). The ambiguity code R has bits 0 and 2 set, corresponding to ambiguity between A and G.
    //           Thus, a bitwise AND with either base's binarized form will be nonzero.
    int match_one(std::vector<unsigned char> const& mask, std::vector<unsigned char> const& refmask) const;
public:
    // Constructor for the BarcodeMatcher
    //  Args:
    //    filename: String containing the path to a tab-delimited file containing at least two columns. The first column is the sequence of each well's barcode, and the second is the name of that well.
    //    max_mismatch: A query barcode will match a reference barcode with at most this many substitutions.
    //    max_n: A query barcode with only be tested for match to a reference barcode if it has at most this many N bases.
    BarcodeMatcher(const std::string& filename, const int &max_mismatch, const int &max_n);

    // Search for the first reference barcode closely matching the query sequence
    //  Args:
    //    seq: The query DNA sequence
    //  Returns:
    //    A pair<int, int> containing the index of the matched barcode and the number of substitutions between it and the query.
    std::pair<int, int> const& match_barcode(std::string const& seq);

    // Access the number of barcodes
    const std::size_t& get_numbc() const {return numbc;}

    // Access the length of each barcode. They are all assumed to be equal in length.
    const std::size_t& get_bclen() const {return bclen;}

    // Access the DNA sequences of the barcodes
    const std::vector<std::string>& get_bcseqs() const {return bcseqs;}

    // Access the names of the wells
    const std::vector<std::string>& get_bcnames() const {return bcnames;}

    // Access the placeholder sequence for noBC (used to fool STAR)
    const std::string& get_nobcseq() const {return nobcseq;}
};

#endif //DRUGSEQ_BARCODEMATCHER_H
