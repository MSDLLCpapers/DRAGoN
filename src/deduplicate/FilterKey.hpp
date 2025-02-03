/*
 * Written by Scott Norton
 * 15-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#ifndef DRUGSEQ_DEDUPLICATE_FILTERKEY_H
#define DRUGSEQ_DEDUPLICATE_FILTERKEY_H

#include <unordered_map>
#include <set>
#include <vector>
#include <bamtools/api/BamAlignment.h>
#include <bamtools/api/SamHeader.h>
#include <drugseq.hpp>

using ChromKey = std::pair<std::string, bool>;

template<> struct std::hash<ChromKey> {
    unsigned long long operator()(ChromKey const& pair) const {
        unsigned long long result = std::hash<std::string>{}(pair.first);
        return (result << 1) ^ pair.second;
    }
};

struct GenomicPos {
    std::string r_id;
    unsigned long long r_pos;
    unsigned long long r_end;
    bool r_strand;

    GenomicPos() = default;
    GenomicPos(const BamTools::RefVector& bam_refs, const BamTools::BamAlignment& bam1) :
            r_id(bam_refs[bam1.RefID].RefName),
            r_pos(bam1.Position),
            r_end(bam1.GetEndPosition()),
            r_strand(bam1.IsReverseStrand())
    {}
    operator std::string() const {
        return r_id + (r_strand ? '-' : '+') + std::to_string(r_pos) + '-' + std::to_string(r_end);
    }
    operator ChromKey() const {
        return std::make_pair(r_id, r_strand);
    }
};

template<> struct std::hash<GenomicPos> {
    unsigned long long operator()(GenomicPos const& pos) const {
        return (std::hash<std::string>{}(pos.r_id) << 1) ^ (std::hash<unsigned long long>{}(pos.r_pos) << 2) ^ (pos.r_strand);
    }
};

struct FilterKey {
    using coord_pair_t = std::pair<const ChromKey, std::vector<unsigned long long>>;
    std::string cb;
    std::string q_name;
    std::unordered_map<ChromKey, std::vector<unsigned long long>> coords;
    unsigned char mapq;
    std::string umi;
    std::set<std::string> genes;
    std::vector<GenomicPos> mx;
    GenomicPos& add_mx(const BamTools::RefVector& bam_refs, const BamTools::BamAlignment& bam1) {
        GenomicPos& mx1 = mx.emplace_back(bam_refs, bam1);
        coords[static_cast<ChromKey>(mx1)].push_back(mx1.r_pos);
        return mx1;
    }
    void reset() {
        cb.clear();
        q_name.clear();
        mapq = 0;
        genes.clear();
        coords.clear();
        mx.clear();
    }
};

#endif //DRUGSEQ_DEDUPLICATE_FILTERKEY_H
