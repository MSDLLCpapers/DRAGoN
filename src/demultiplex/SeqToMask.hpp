#ifndef DRUGSEQ_DEMULTIPLEX_SEQTOMASK_H
#define DRUGSEQ_DEMULTIPLEX_SEQTOMASK_H

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#include <vector>
#include <algorithm>
#include <bamtools/api/BamConstants.h>

static inline unsigned char base_to_mask(char c) {
    switch (c) {
    case (BamTools::Constants::BAM_DNA_EQUAL):
        return BamTools::Constants::BAM_BASECODE_EQUAL;
    case (BamTools::Constants::BAM_DNA_A):
        return BamTools::Constants::BAM_BASECODE_A;
    case (BamTools::Constants::BAM_DNA_C):
        return BamTools::Constants::BAM_BASECODE_C;
    case (BamTools::Constants::BAM_DNA_M):
        return BamTools::Constants::BAM_BASECODE_M;
    case (BamTools::Constants::BAM_DNA_G):
        return BamTools::Constants::BAM_BASECODE_G;
    case (BamTools::Constants::BAM_DNA_R):
        return BamTools::Constants::BAM_BASECODE_R;
    case (BamTools::Constants::BAM_DNA_S):
        return BamTools::Constants::BAM_BASECODE_S;
    case (BamTools::Constants::BAM_DNA_V):
        return BamTools::Constants::BAM_BASECODE_V;
    case (BamTools::Constants::BAM_DNA_T):
        return BamTools::Constants::BAM_BASECODE_T;
    case (BamTools::Constants::BAM_DNA_W):
        return BamTools::Constants::BAM_BASECODE_W;
    case (BamTools::Constants::BAM_DNA_Y):
        return BamTools::Constants::BAM_BASECODE_Y;
    case (BamTools::Constants::BAM_DNA_H):
        return BamTools::Constants::BAM_BASECODE_H;
    case (BamTools::Constants::BAM_DNA_K):
        return BamTools::Constants::BAM_BASECODE_K;
    case (BamTools::Constants::BAM_DNA_D):
        return BamTools::Constants::BAM_BASECODE_D;
    case (BamTools::Constants::BAM_DNA_B):
        return BamTools::Constants::BAM_BASECODE_B;
    case (BamTools::Constants::BAM_DNA_N):
    default:
        return BamTools::Constants::BAM_BASECODE_N;
    }
}

// Binarize DNA string. The output is populated with 4-bit characters.
//  Args:
//    seq: DNA sequence, composed of base or ambiguity codes.
//    vec: Vector of unsigned char constituting the output
static inline void seq_to_mask(const std::string& seq, std::vector<unsigned char>& vec) {
    std::ranges::transform(seq, vec.begin(), base_to_mask);
}

#endif //DRUGSEQ_DEMULTIPLEX_SEQTOMASK_H
