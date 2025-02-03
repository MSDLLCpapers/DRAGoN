/*
 * Written by Scott Norton
 * 15-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#ifndef DRUGSEQ_DEDUPLICATE_DEDUPLICATEWORKER_H
#define DRUGSEQ_DEDUPLICATE_DEDUPLICATEWORKER_H

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <filesystem>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/SparseExtra>
#include "SamReader.hpp"
#include <bamtools/api/BamWriter.h>
#include <bamtools/api/BamAlignment.h>
#include "DuplicateDetector.hpp"
#include "MultimapHandlers.hpp"

namespace fs = std::filesystem;

template<typename Scalar> using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::ColMajor, long long>;

typedef SparseMatrix<unsigned long long> unique_mtx;
typedef SparseMatrix<double> multimaps_mtx;

class DeduplicateWorker {
    SamReader reader; // Handles read operations from the source BAM
    BamTools::BamWriter writer; // Handles write operations to the dest BAM
    DuplicateDetector detector; // Steam FilterKeys through this to detect duplicates
    std::vector<std::vector<unsigned long long>> stats; // Records QC and UMI statistics
    std::unordered_map<std::string, std::pair<int, std::string>> barcodes; // Read from barcodes.txt
    std::unordered_map<std::string, std::pair<unsigned long long, std::string>> features; // Read from GTF
    std::map<std::set<std::string>, std::size_t> multimapped_reads; // Defer counting reads mapping to two or more genes
    std::vector<std::pair<std::vector<unsigned long long>, std::size_t>> mmaps_flat; // std::flat_map does not exist in C++17
    unique_mtx unique_counts; // Integer sparse matrix counting reads mapping to a single gene
    std::unordered_map<MultimapStrategy, multimaps_mtx> matrices;
    std::size_t numbc, numft;
    const int em_iter_max;
    const double em_max_diff;
    const double em_small_thresh;
    const bool keep_bam;

    void read_barcodes(const std::string& barcodes_in);
    void read_features(const std::string& gtf_in);
    void preallocate_matrices(const std::unordered_set<MultimapStrategy>& mmap_strategies);
    void handle_multimapped_reads(const int& cbidx);
    void handle_multimaps_Uniform(const int& cbidx, multimaps_mtx& mtx);
    void handle_multimaps_PropUnique(const int& cbidx, multimaps_mtx& mtx);
    void handle_multimaps_EM(const int& cbidx, multimaps_mtx& mtx);
    void handle_multimaps_Rescue(const int& cbidx, multimaps_mtx& mtx);
    void dump_barcodes(const fs::path& fname) const;
    void dump_features(const fs::path& fname) const;
    void dump_stats(const fs::path& fname) const;

    template<typename Vt = double>
    static void dump_matrix(const fs::path& fname, const SparseMatrix<Vt>& mtx) {
        Eigen::saveMarket(mtx, fname);
    }

public:
    DeduplicateWorker(
        const std::vector<std::string>& cli,
        const std::string& bamfilename,
        const unsigned long long& threads,
        const std::string& barcodes_in,
        const std::string& gtf_in,
        const std::unordered_set<MultimapStrategy>& mmap_strategies,
        const unsigned long long& umilen = 10,
        const unsigned long long& max_distance = 2,
        const unsigned long long& max_lookback = 10,
        const unsigned long long& em_iter_max = 100,
        const double& em_max_diff = 0.01,
        const double& em_small_thresh = 0.01,
        const bool& keep_bam = false
    );
    int run();
    void save_output(fs::path const& outdir) const;
};

#endif //DRUGSEQ_DEDUPLICATE_DEDUPLICATEWORKER_H
