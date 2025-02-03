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
#include <array>
#include <algorithm>
#include "ReadPair.hpp"
#include "QCStats.hpp"
#include "QualityFilter.hpp"
#include "SeqToMask.hpp"

QCFiltrator::QCFiltrator(const int &min_len, const int &min_qual, const int &phred, const int &umi_qual_bases, const int &homop_len, const double &max_AT, const double &max_GC, const std::string &adapter) :
    min_len(min_len),
    min_qual(min_qual + phred),
    umi_qual_bases(umi_qual_bases),
    homop_len(homop_len),
    max_AT(max_AT),
    max_GC(max_GC),
    adapter(adapter)
{
    enable_length = true;
    if (min_qual > 0) {
        enable_trimming = true;
        if (umi_qual_bases > 0) {
            enable_umiqual = true;
        }
    }
    if (max_AT < 1.0 && max_GC < 1.0) {
        enable_homogeny = true;
    }
    if (adapter.length() >= 5) {
        enable_adapter_trimming = true;
        adapter_alignment_matrix.resize(adapter.length() + 1);
    }
    if (homop_len > 5) {
        enable_homopolymer = true;
        homopolymers.at(0).assign(homop_len, 'A');
        homopolymers.at(1).assign(homop_len, 'C');
        homopolymers.at(2).assign(homop_len, 'G');
        homopolymers.at(3).assign(homop_len, 'T');
        homopolymers.at(4).assign(homop_len, 'N');
    }

    enable_duplicate_filter = false;
}

std::optional<CutAdaptDPMatrixCell> QCFiltrator::cutadapt_match_adapter(ReadPair const& rp) {
    // C++ implementation of CutAdapt's trimming algorithm
    // reference = adapter
    // query = rp.cdna_seq

    // Hardcode as a 3' matcher and trimmer
    // BackAdapter flags: QUERY_START | QUERY_STOP | REFERENCE_END
    int i, j;
    int k = 1;  // max mismatches
    int max_n = rp.cdna_seq.length(), max_m = adapter.length();
    int last = std::min(max_m, k + 1);
    int last_filled_i;
    // int min_n = 0;
    CutAdaptDPMatrixCell best {.cost = static_cast<int>(max_n + max_m + 1)};
    int best_ref_stop = max_m;

    // Before first base of adapter
    std::iota(adapter_alignment_matrix.begin(), adapter_alignment_matrix.end(), CutAdaptDPMatrixCell{0, 0, 0});
    // Nth base
    for (j = 1; j <= max_n; ++j) {
        CutAdaptDPMatrixCell diag_entry = adapter_alignment_matrix.front();  // because we only keep one row in memory, save the entry to the upper left
        ++adapter_alignment_matrix.front().origin;
        for (i = 1; i <= last; ++i) {
            if ((base_to_mask(adapter.at(i - 1)) & base_to_mask(rp.cdna_seq.at(j - 1))) != 0) {  // IUPAC match
                ++diag_entry.score;  // match score increment
            } else {
                std::array<int, 3> costs {
                    diag_entry.cost + 1,                             // mismatch
                    adapter_alignment_matrix.at(i).cost + 2,      // deletion
                    adapter_alignment_matrix.at(i - 1).cost + 2,  // insertion
                };
                auto best = std::min_element(costs.begin(), costs.end());
                diag_entry.cost = *best;
                --diag_entry.score;  // fixed penalty for mismatch, insertion, and deletion are all 1
                switch (best - costs.begin()) {
                case 0:  // mismatch
                    break;
                case 1:  // deletion
                    diag_entry.origin = adapter_alignment_matrix.at(i).origin;  // directly above
                    break;
                case 2:  // insertion
                    diag_entry.origin = adapter_alignment_matrix.at(i - 1).origin;  // directly to the left
                    break;
                }
            }
            std::swap(diag_entry, adapter_alignment_matrix.at(i));
        }
        last_filled_i = last;
        // restruct last to be at most k
        // could result in last = -1 but we increment directly after
        while (last >= 0 && adapter_alignment_matrix.at(last).cost > k) {--last;}
        if (last < max_m) {
            ++last;
        } else {  // if (stop_in_query)
            const CutAdaptDPMatrixCell& curr = adapter_alignment_matrix.back();
            int length = i + std::min(0, curr.origin);
            int best_length = best_ref_stop + std::min(best.origin, 0);
            if ((length >= 5 && curr.cost <= k) && (
                best.cost == max_n + max_m + 1 ||
                (curr.origin <= best.origin + max_m / 2 && curr.score >= best.score) ||
                (length >= best_length && curr.score >= best.score)
            )) {
                best = curr;
                if (curr.cost == 0 && curr.origin >= 0) {
                    break;
                }
            }
        }
    }
    // if (max_n == n)
    int first_i = 0;  // if stop_in_reference

    // Find the best starting point within the cDNA
    for (i = last_filled_i; i >= first_i; --i) {
        const CutAdaptDPMatrixCell& curr = adapter_alignment_matrix.at(i);
        int length = i + std::min(0, curr.origin);
        int best_length = best_ref_stop + std::min(best.origin, 0);
        if ((length >= 5 && curr.cost <= k) && (
            best.cost == max_n + max_m + 1 ||
            (curr.origin <= best.origin + max_m / 2 && curr.score >= best.score) ||
            (length >= best_length && curr.score >= best.score)
        )) {
            best = curr;
            best_ref_stop = i;
        }
    }

    if (best.cost == max_n + max_m + 1) {
        return {};
    }
    return best;
}

QCFailReason QCFiltrator::operator()(ReadPair& rp) {
    if (enable_umiqual) {
        if (std::count_if(rp.umiqual.cbegin(), rp.umiqual.cend(), [&](const char &c) {
            return c < min_qual;
        }) > umi_qual_bases) {
            return QCFailReason::FAIL_UMI_QUALITY;
        }
    }
    std::size_t cdna_length = rp.cdna_seq.length();
    if (enable_homogeny) {
        int count_AT = 0, count_GC = 0;
        for (const char &c : rp.cdna_seq) {
            switch (c) {
            case 'A':
            case 'T':
                ++count_AT;
                break;
            case 'C':
            case 'G':
                ++count_GC;
                break;
            }
        }
        if (count_AT >= max_AT * cdna_length) {
            return QCFailReason::FAIL_AT_CONTENT;
        }
        if (count_GC >= max_GC * cdna_length) {
            return QCFailReason::FAIL_GC_CONTENT;
        }
    }
    if (enable_length && cdna_length < min_len) {
        return QCFailReason::FAIL_LENGTH;
    }
    if (enable_homopolymer) {
        std::size_t pos;
        for (const std::string& homop : homopolymers) {
            pos = rp.cdna_seq.find(homop);
            if (pos != std::string::npos) {
                if (pos < min_len) {
                    return QCFailReason::FAIL_HOMOPOLYMER;
                }
                rp.trim_3prime(cdna_length = pos);
            }
        }
    }
    if (enable_adapter_trimming) {
        std::optional<CutAdaptDPMatrixCell> best = cutadapt_match_adapter(rp);
        if (best) {
            // Found a match. If removing it won't result in a too-short read, lop it off
            if (best->origin < min_len) {
                return QCFailReason::FAIL_ADAPTER_TRIMMING;
            }
            rp.trim_3prime(cdna_length = best->origin);
        }
    }
    if (enable_trimming) {
        // C++ implementation of BWA Nextseq trimming algorithm
        int s = 0;
        int max_qual = 0;
        int max_i = cdna_length;
        for (int i = cdna_length - 1; i >= 0; --i) {
            char q = rp.cdna_qual.at(i);
            if (rp.cdna_seq.at(i) == 'G') {
                q = min_qual - 1;
            }
            s += min_qual - q;
            if (s < 0) {
                break;
            }
            if (s > max_qual) {
                max_qual = s;
                max_i = i;
            }
        }
        if (max_i < min_len) {
            return QCFailReason::FAIL_QUALITY_TRIMMING;
        }
        rp.trim_3prime(max_i);
    }
    // Rudimentary duplicate detection
    if (enable_duplicate_filter) {
        if (seen.find(rp) != seen.end()) {
            return QCFailReason::FAIL_DUPLICATE;
        }
        seen.insert(rp);
    }

    return QCFailReason::PASS;
}
