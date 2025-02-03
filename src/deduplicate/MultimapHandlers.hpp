/*
 * Written by Scott Norton
 * 15-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#ifndef DRUGSEQ_DEDUPLICATE_MULTIMAPHANDLERS_H
#define DRUGSEQ_DEDUPLICATE_MULTIMAPHANDLERS_H

enum MultimapStrategy {
    UNIQUE = 0,
    UNIFORM,
    PROPUNIQUE,
    EM,
    RESCUE
};

#endif //DRUGSEQ_DEDUPLICATE_MULTIMAPHANDLERS_H
