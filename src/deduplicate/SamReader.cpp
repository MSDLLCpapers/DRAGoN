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
#include <iterator>
#include <iostream>
#include <algorithm>
#include <deque>
#include <drugseq.hpp>
#include "DedupStats.hpp"
#include <bamtools/api/BamAlignment.h>
#include <bamtools/api/BamReader.h>
#include <bamtools/api/SamHeader.h>
#include "SamReader.hpp"

using namespace std::chrono_literals;

SamReader::SamReader(const std::string& filename, const unsigned long long threads) {
    Open(filename);
}


DedupStat SamReader::get_read(FilterKey& out, BamTools::BamAlignment& primary) {
    BamTools::BamReader& parent = *this;
    return get_read_impl(parent, out, primary);
}

DedupStat SamReader::get_read_impl(BamTools::BamReader& samfile, FilterKey& out, BamTools::BamAlignment& primary) {
    static BamTools::BamAlignment bam1;
    static bool is_init = false;
    static bool is_eof = false;
    static bool is_error = false;
    std::deque<std::string> mx, xt, xs;
    int bamread_result;

    out.reset();

    if (!is_init) {
        if (!samfile.GetNextAlignment(bam1)) {
            if (samfile.GetErrorString().empty()) {
                is_eof = true;
                return DedupStat::FAIL_EOF;
            } else {
                is_error = true;
                return DedupStat::FAIL_HTS_ERROR;
            }
        }
        is_init = true;
    } else if (is_eof) {
        return DedupStat::FAIL_EOF;
    } else if (is_error) {
        return DedupStat::FAIL_HTS_ERROR;
    }
    out.q_name = bam1.Name;
    if (!bam1.GetTag<std::string>("CB", out.cb)) {
        is_error = true;
        return DedupStat::FAIL_HTS_ERROR;
    }
    if (!bam1.GetTag<std::string>("UR", out.umi)) {
        is_error = true;
        return DedupStat::FAIL_HTS_ERROR;
    }
    unsigned char nh;
    if (!bam1.GetTag<unsigned char>("NH", nh)) {
        is_error = true;
        return DedupStat::FAIL_HTS_ERROR;
    }
    for (int i = 0; i < nh; ++i) {
        std::string cur_xt {"-"}, cur_xs, cur_mx;
        if (bam1.GetTag<std::string>("XT", cur_xt)) {
            drugseq::strsplit(std::inserter(out.genes, out.genes.end()), cur_xt, ",");
        }
        bam1.GetTag<std::string>("XS", cur_xs);
        cur_mx = static_cast<std::string>(out.add_mx(samfile.GetReferenceData(), bam1));
        if (!bam1.IsPrimaryAlignment()) {
            xt.push_back(cur_xt);
            xs.push_back(cur_xs);
            mx.push_back(cur_mx);
        } else {
            primary = bam1;
            xt.push_front(cur_xt);
            xs.push_front(cur_xs);
            //mx.push_front(cur_mx);
        }
        if (!samfile.GetNextAlignment(bam1)) {
            if (samfile.GetErrorString().empty()) {
                is_eof = true;
                break;
            } else {
                is_error = true;
                return DedupStat::FAIL_HTS_ERROR;
            }
        } else if (i != nh - 1 && bam1.Name != out.q_name) {
            is_error = true;
            std::cerr << "FATAL: BAM file not collated\n";
            return DedupStat::FAIL_HTS_ERROR;
        }
    }
    if (nh == 0) {
        primary = bam1;
        if (!samfile.GetNextAlignment(bam1)) {
            if (samfile.GetErrorString().empty()) {
                is_eof = true;
            } else {
                is_error = true;
                return DedupStat::FAIL_HTS_ERROR;
            }
        }
    }
    if (!primary.IsPrimaryAlignment()) {
        is_error = true;
        std::cerr << "FATAL: No primary alignment in read group\n";
        return DedupStat::FAIL_HTS_ERROR;
    }
    if (nh > 1) {
        if (
            !primary.EditTag("MX", "Z", drugseq::strjoin(mx, ",")) ||
            (!out.genes.empty() && !primary.EditTag("XN", "I", static_cast<uint32_t>(out.genes.size()))) ||
            !primary.EditTag("XS", "Z", drugseq::strjoin(xs, ",")) ||
            !primary.EditTag("XT", "Z", drugseq::strjoin(xt, ","))
        ) {
            return DedupStat::FAIL_HTS_ERROR;
        }
    }

    // QC filtering
    std::string qcfail_reason;
    return
        (
            primary.IsFailedQC()
            || !primary.GetTag<std::string>("QC", qcfail_reason)
            || qcfail_reason != "PASS"
        )
            ? DedupStat::DISCARD_QCFAIL
            : !primary.IsMapped()
                ? DedupStat::DISCARD_UNMAPPED
                : std::find(xs.begin(), xs.end(),"Assigned") == xs.end()
                ? DedupStat::DISCARD_UNASSIGNED
                : DedupStat::PASSED_QC;
}
