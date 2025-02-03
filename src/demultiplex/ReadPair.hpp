/*
 * Written by Scott Norton
 * 12-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#ifndef DRUGSEQ_READPAIR_H
#define DRUGSEQ_READPAIR_H

#include <string>
#include <stdexcept>

struct ReadPair {
    std::string name;
    std::string cdna_seq;
    std::string cdna_qual;
    std::string bcseq;
    std::string umiseq;
    std::string bcqual;
    std::string umiqual;

    ReadPair& trim_3prime(size_t pos) {
        if (pos != std::string::npos && pos > cdna_seq.length()) {
            throw std::range_error("ReadPair::trim_3prime");
        }
        cdna_seq = cdna_seq.substr(0, pos);
        cdna_qual = cdna_qual.substr(0, pos);
        return *this;
    }

    bool operator==(const ReadPair& rhs) const {
        return cdna_seq == rhs.cdna_seq && bcseq == rhs.bcseq && umiseq == rhs.umiseq;
    }
};

template<> struct std::hash<ReadPair> {
    std::size_t operator()(const ReadPair& rp) const {
        std::hash<std::string> strhash{};
        return (strhash(rp.cdna_seq) << 6) ^ (strhash(rp.bcseq) << 1) ^ strhash(rp.umiseq);
    }
};

#endif //DRUGSEQ_READPAIR_H
