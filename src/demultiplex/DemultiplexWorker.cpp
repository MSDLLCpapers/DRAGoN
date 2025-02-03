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
#include <vector>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <bamtools/api/BamWriter.h>
#include <bamtools/api/SamHeader.h>
#include <bamtools/api/SamReadGroup.h>
#include <bamtools/api/BamAlignment.h>
#include <drugseq.hpp>
#include "QCStats.hpp"
#include "BarcodeMatcher.hpp"
#include "QualityFilter.hpp"
#include "FastqReader.hpp"
#include "DemultiplexWorker.hpp"

using namespace std::literals;

template <typename Stream, typename TimePoint, typename Clock = typename TimePoint::clock>
Stream& operator<<(Stream& strm, const TimePoint& timepoint) {
    auto ctime = Clock::to_time_t(timepoint);
    strm << std::put_time(std::localtime(&ctime), "%F %X");
    return strm;
}

MainWorker::MainWorker(
    const std::string& barcode_filename,
    const std::string& bamfpref,
    const std::vector<std::string>& cli,
    const int& umilen,
    const int& split,
    const int& max_mismatch,
    const int& max_n,
    const int& min_len,
    const int& min_qual,
    const int& phred,
    const int& umi_qual_bases,
    const int& homop_len,
    const double& max_AT,
    const double& max_GC,
    const int& bam_threads,
    const std::string& adapter,
    const bool write_qcfail
) :
    split(split),
    bc(barcode_filename, max_mismatch, max_n),
    qc(min_len, min_qual, phred, umi_qual_bases, homop_len, max_AT, max_GC, adapter),
    stat(bc, split ? split : bc.get_numbc(), max_mismatch),
    reader(bc.get_bclen(), umilen),
    write_qcfail(write_qcfail)
{
    BamTools::SamHeader header {};
    header.Version = "1.6";
    header.SortOrder = "unsorted";
    BamTools::SamProgram program {cli.front()};
    program.Name = cli.front();
    program.CommandLine = drugseq::shlexjoin(cli);
    program.Version = DRUGSEQ_VERSION_STR;
    header.Programs.Add(program);
    std::string bcidx_s {"0000"};
    for (int i = 0; i < stat.max_file; ++i) {
        header.ReadGroups[bcidx_s].SequencingTechnology = "ILLUMINA";
        header.ReadGroups[bcidx_s].Sample = bamfpref;
        if (split == 0) {
            header.ReadGroups[bcidx_s].Description = bc.get_bcseqs().at(i);
        } else {
            int bcidx_lo = i * bc.get_numbc() / split;
            int bcidx_hi = (i + 1) * bc.get_numbc() / split;
            header.ReadGroups[bcidx_s].Description = drugseq::strjoin(bc.get_bcseqs().begin() + bcidx_lo, bc.get_bcseqs().begin() + bcidx_hi, "-"s);
        }
        for (auto it = bcidx_s.rbegin(); it != bcidx_s.rend(); ++it) {
            char& digit = *it;
            if (++digit <= '9') break;
            digit = '0';
        }
    }
    header.ReadGroups["noBC"].SequencingTechnology = "ILLUMINA";
    header.ReadGroups["noBC"].Sample = bamfpref;
    header.ReadGroups["noBC"].Description = std::string(bc.get_bclen(), 'N');
    std::string header_s = header.ToString();

    #if BAM_WRITE_MULTIPLE
    outbam.reserve(stat.max_file + 1);
    bcidx_s = "0000";
    for (int i = 0; i < stat.max_file; ++i) {
        outbam.emplace_back().Open(bamfpref + "_" + bcidx_s + ".bam", header_s, {});
        for (auto it = bcidx_s.rbegin(); it != bcidx_s.rend(); ++it) {
            char& digit = *it;
            if (++digit <= '9') break;
            digit = '0';
        }
    }
    outbam.emplace_back().Open(bamfpref + "_noBC.bam", header_s, {});
    #else
    outbam.reserve(2);
    outbam.emplace_back().Open(bamfpref + ".bam", header_s, {});
    #endif //BAM_WRITE_MULTIPLE
    outbam.emplace_back().Open(bamfpref + "_qcFail.bam", header_s, {});
}

bool MainWorker::write_bam(const ReadPair& rp, const std::string& bcseq, const int fileidx, const int mismatch, const QCFailReason qcfail) {
    static BamTools::BamAlignment bam1;
    std::string bi_tag = "noBC";
    if (fileidx >= 0) {
        std::stringstream ss;
        ss << std::setw(4) << std::setfill('0') << fileidx;
        bi_tag = ss.str();
    }
    bam1.Name = rp.name;
    bam1.QueryBases = rp.cdna_seq;
    bam1.Qualities = rp.cdna_qual;
    bam1.SetIsMapped(false);
    bam1.SetIsFailedQC(qcfail != QCFailReason::PASS);
    auto& bamfile = qcfail == QCFailReason::PASS ?
        #if BAM_WRITE_MULTIPLE
        outbam.at((fileidx == -1 ? stat.max_file : fileidx))
        #else
        outbam.front()
        #endif //BAM_WRITE_MULTIPLE
    : outbam.back();
    return bam1.EditTag("CR", "Z", rp.bcseq) &&
           bam1.EditTag("UR", "Z", rp.umiseq) &&
           bam1.EditTag("CY", "Z", rp.bcqual) &&
           bam1.EditTag("UY", "Z", rp.umiqual) &&
           bam1.EditTag("CB", "Z", bcseq) &&
           bam1.EditTag("UB", "Z", rp.umiseq) &&
           bam1.EditTag("BI", "Z", bi_tag) &&
           bam1.EditTag("BM", "i", mismatch) &&
           bam1.EditTag("QC", "Z", QCFailNames[qcfail]) &&
           bam1.EditTag("RG", "Z", bi_tag) &&
           bamfile.SaveAlignment(bam1);
}

int MainWorker::run(const std::string& r1fq, const std::string& r2fq) {
    if (!reader.open(r1fq, r2fq)) {
        throw std::runtime_error("Unable to open FASTQ file pair for reading");
    }
    ReadPair rp;
    std::cerr << std::chrono::system_clock::now() << " [INFO] Begin\n";
    int total_in = 0;
    while (reader.next_read_pair(rp)) {
        int bc_idx;
        std::pair<int, int> const& match_result = bc.match_barcode(rp.bcseq);
        QCFailReason qcfail = qc(rp);
        int fileidx = match_result.first;
        if (match_result.first >= 0 && split != 0) {
            fileidx = fileidx * split / bc.get_numbc();
        }
        stat.record(rp, qcfail, match_result.first, match_result.second, fileidx);
        if ((++total_in % 1000000) == 0) {
            std::cerr << std::chrono::system_clock::now() << " [INFO] Processed " << total_in << " reads...\n";
        }
        if (qcfail == QCFailReason::PASS || write_qcfail) {
            if (!write_bam(rp, match_result.first >= 0 ? bc.get_bcseqs().at(match_result.first) : bc.get_nobcseq(), fileidx, match_result.second, qcfail)) {
                throw std::runtime_error("SAM write error");
            }
        }
    }
    std::cerr << std::chrono::system_clock::now() << " [INFO] Finished! " << total_in << " reads processed\n";
    return 0;
}

void MainWorker::dump_stats(std::ostream& strm) {
    strm << stat;
}
