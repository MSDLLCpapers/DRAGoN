/*
 * Written by Scott Norton
 * 12-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#ifndef DRUGSEQ_FASTQREADER_H
#define DRUGSEQ_FASTQREADER_H

#include <string>
#include <kseq++/seqio.hpp>
#include "ReadPair.hpp"

class FastqReader {
    int bclen, umilen, phred = 33;
    std::unique_ptr<klibpp::SeqStreamIn> fastq1, fastq2;

public:
    FastqReader(const int &bclen, const int &umilen);
    FastqReader(const std::string &r1fqname, const std::string &r2fqname, const int &bclen, const int &umilen);
    bool open(const std::string &r1fqname, const std::string &r2fqname);
    bool next_read_pair(ReadPair& rp);
};

#endif //DRUGSEQ_FASTQREADER_H
