/*
 * Written by Scott Norton
 * 15-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#ifndef DRUGSEQ_DEDUPLICATE_SAMREADER_H
#define DRUGSEQ_DEDUPLICATE_SAMREADER_H

#ifndef SAMREADER_THREAD_IMPL
#define SAMREADER_THREAD_IMPL 0
#endif

#include <string>
#include "FilterKey.hpp"
#include "DedupStats.hpp"
#include <drugseq.hpp>
#include <bamtools/api/BamAlignment.h>
#include <bamtools/api/BamReader.h>
#include <bamtools/api/SamHeader.h>

class SamReader : public BamTools::BamReader
{
    static DedupStat get_read_impl(BamTools::BamReader& samfile, FilterKey& out, BamTools::BamAlignment& primary);
public:
    SamReader(const std::string& filename, const unsigned long long threads);
    DedupStat get_read(FilterKey& out, BamTools::BamAlignment& primary);
};

#endif //DRUGSEQ_DEDUPLICATE_SAMREADER_H
