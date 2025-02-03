/*
 * Written by Scott Norton
 * 12-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#ifndef DRUGSEQ_DEMULTIPLEXWORKER_H
#define DRUGSEQ_DEMULTIPLEXWORKER_H

#define BAM_WRITE_MULTIPLE        (0)

#include <string>
#include <vector>
#include <bamtools/api/BamWriter.h>
#include "QCStats.hpp"
#include "BarcodeMatcher.hpp"
#include "QualityFilter.hpp"
#include "FastqReader.hpp"

class MainWorker {
    const int& split;
    BarcodeMatcher bc;
    QCFiltrator qc;
    Statistics stat;
    FastqReader reader;
    std::vector<BamTools::BamWriter> outbam;
    bool write_qcfail = true;
public:
    MainWorker(
        const std::string& barcode_filename,
        const std::string& bamfname,
        const std::vector<std::string>& cli,
        const int& umilen = 10,
        const int& split = 0,
        const int& max_mismatch = 1,
        const int& max_n = 1,
        const int& min_len = 20,
        const int& min_qual = 20,
        const int& phred = 33,
        const int& umi_qual_bases = 6,
        const int& homop_len = 10,
        const double& max_AT = 0.9,
        const double& max_GC = 0.9,
        const int& bamthreads = 1,
        const std::string& adapter = "",
        const bool write_qcfail = true
    );
    int run(const std::string& r1fq, const std::string& r2fq);
    void dump_stats(std::ostream& strm);
private:
    bool write_bam(const ReadPair& rp, const std::string& bcseq, const int fileidx, const int mismatch, const QCFailReason qcfail);
};

#endif //DRUGSEQ_DEMULTIPLEXWORKER_H
