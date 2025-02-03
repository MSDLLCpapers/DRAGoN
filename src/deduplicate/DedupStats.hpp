/*
 * Written by Scott Norton
 * 15-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#ifndef DRUGSEQ_DEDUPLICATE_DEDUPSTATS_H
#define DRUGSEQ_DEDUPLICATE_DEDUPSTATS_H

enum DedupStat {
    TOTAL_IN = 0,
    DISCARD_QCFAIL,
    DISCARD_UNMAPPED,
    DISCARD_UNASSIGNED,
    PASSED_QC,
    N_MULTIMAPPED,
    UMI_COUNT,
    UNIQUE_COUNT,
    UNIQUE_UMIS,
    MAX,

    FAIL_EOF = -1,
    FAIL_HTS_ERROR = -2
};

#endif //DRUGSEQ_DEDUPLICATE_DEDUPSTATS_H
