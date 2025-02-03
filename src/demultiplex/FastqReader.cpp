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
#include <drugseq.hpp>
#include <kseq++/seqio.hpp>
#include "ReadPair.hpp"
#include "FastqReader.hpp"

FastqReader::FastqReader(const int &bclen, const int &umilen) : bclen(bclen), umilen(umilen)
{}

FastqReader::FastqReader(const std::string &r1fqname, const std::string &r2fqname, const int &bclen, const int &umilen) : bclen(bclen), umilen(umilen)
{
    open(r1fqname, r2fqname);
}

bool FastqReader::open(const std::string &r1fqname, const std::string &r2fqname) {
    fastq1 = std::make_unique<klibpp::SeqStreamIn>(r1fqname.c_str());
    fastq2 = std::make_unique<klibpp::SeqStreamIn>(r2fqname.c_str());
    return !!*fastq1 && !!*fastq2;
}

bool FastqReader::next_read_pair(ReadPair& rp) {
    klibpp::KSeq r1, r2;
    (*fastq1) >> r1;
    (*fastq2) >> r2;
    if (r1.name.empty() || r2.name.empty()) {
        return false;
    }
    if (r1.seq.length() < bclen + umilen) {
        throw std::runtime_error("R1: Invalid fastq record: sequence too short (" + std::to_string(r1.seq.length()) + " < " + std::to_string(bclen + umilen) + ")");
    }
    rp.name = r1.name;
    rp.cdna_seq = r2.seq;
    rp.cdna_qual = r2.qual;
    rp.bcseq = r1.seq.substr(0, bclen);
    rp.umiseq = r1.seq.substr(bclen, umilen);
    rp.bcqual = r1.qual.substr(0, bclen);
    rp.umiqual = r1.qual.substr(bclen, umilen);
    return true;
}
