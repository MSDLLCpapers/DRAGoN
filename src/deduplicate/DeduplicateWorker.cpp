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
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <numeric>
#include <filesystem>
#include <iostream>
#include <optional>
#include <utility>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/SparseExtra>
#include "SamReader.hpp"
#include "DuplicateDetector.hpp"
#include "MultimapHandlers.hpp"
#include "DeduplicateWorker.hpp"
#include <drugseq.hpp>

namespace fs = std::filesystem;
using namespace std::chrono_literals;

DeduplicateWorker::DeduplicateWorker(
    const std::vector<std::string>& cli,
    const std::string& bamfilename,
    const unsigned long long& threads,
    const std::string& barcodes_in,
    const std::string& gtf_in,
    const std::unordered_set<MultimapStrategy>& mmap_strategies,
    const unsigned long long& umilen,
    const unsigned long long& max_distance,
    const unsigned long long& max_lookback,
    const unsigned long long& em_iter_max,
    const double& em_max_diff,
    const double& em_small_thresh,
    const bool& keep_bam
) :
        reader(bamfilename, threads),
        detector(umilen, max_distance, max_lookback),
        numbc(0), numft(0),
        em_iter_max(em_iter_max),
        em_max_diff(em_max_diff),
        em_small_thresh(em_small_thresh),
        keep_bam(keep_bam)
{
    if (keep_bam) {
        BamTools::SamHeader header = reader.GetHeader();
        BamTools::SamProgram program (cli[0]);
        program.Name = cli[0];
        program.Version = DRUGSEQ_VERSION_STR;
        program.CommandLine = drugseq::shlexjoin(cli);
        header.Programs.Add(program);
        writer.Open(bamfilename.substr(0, bamfilename.rfind('.')) + "_MarkDups.bam", header, reader.GetReferenceData());
    }
    read_barcodes(barcodes_in);
    read_features(gtf_in);
    preallocate_matrices(mmap_strategies);
}

void DeduplicateWorker::read_barcodes(const std::string& barcodes_in) {
    std::ifstream bcfile(barcodes_in);
    std::string line;
    std::cerr << drugseq::put_time() << ": Reading barcodes " << barcodes_in << std::endl;
    std::size_t bclen;
    while (std::getline(bcfile, line)) {
        std::size_t space1 = line.find('\t'), space2 = line.find('\t', space1 + 1);
        if (space1 == std::string::npos) {
            break;
        }
        bclen = space1;
        barcodes[line.substr(0, space1)] = std::make_pair(numbc++, line.substr(space1 + 1, space2 - space1 - 1));
    }
    if (barcodes.empty()) {
        throw std::runtime_error("no barcodes found");
    }
    barcodes[std::string(bclen, 'N')] = std::make_pair(-1, "noBC");
}

void DeduplicateWorker::read_features(const std::string& gtf_in) {
    std::cerr << drugseq::put_time() << ": Reading features " << gtf_in << std::endl;
    drugseq::Gtf gtf {gtf_in};
    for (drugseq::GtfRecord record : gtf.lazy_load("gene")) {
        const std::string &gene_id = std::get<std::string>(record.attr.at("gene_id"));
        std::optional<std::string> gene_name;
        try {
            gene_name = std::get<std::string>(record.attr.at("gene_name"));
        } catch (std::out_of_range) {}
        features[gene_id] = std::make_pair(numft++, gene_name.value_or(""));
    }
    if (features.empty()) {
        throw std::runtime_error("failed to load GTF features");
    }
    std::cerr << drugseq::put_time() << ": Read " << numft << " features " << gtf_in << std::endl;

}

void DeduplicateWorker::preallocate_matrices(const std::unordered_set<MultimapStrategy>& mmap_strategies) {
    // Resize stats matrix
    stats.resize(numbc + 1);
    for (std::vector<unsigned long long>& row : stats) {
        row.resize(DedupStat::MAX);
    }
    // Set sparse matrix dims
    unique_counts.resize(numft, numbc);
    unique_counts.uncompress();
    for (const MultimapStrategy &strat : mmap_strategies) {
        matrices[strat] = multimaps_mtx(numft, numbc);
        matrices[strat].uncompress();
    }
}

void DeduplicateWorker::handle_multimaps_Uniform(const int& cbidx, multimaps_mtx& mtx) {
    for (unique_mtx::InnerIterator it (unique_counts, cbidx); it; ++it) {
        mtx.coeffRef(it.row(), cbidx) = it.value();
    }
    for (const std::pair<const std::vector<unsigned long long>, std::size_t>& umiset : mmaps_flat) {
        // Expand it out to force floating-point maths
        double factor = 1.0;
        factor *= umiset.second;
        factor /= umiset.first.size();
        for (const unsigned long long &ftidx : umiset.first) {
            mtx.coeffRef(ftidx, cbidx) += factor;
        }
    }
}

void DeduplicateWorker::handle_multimaps_PropUnique(const int& cbidx, multimaps_mtx& mtx) {
    std::vector<double> uniques(numft);
    for (unique_mtx::InnerIterator it (unique_counts, cbidx); it; ++it) {
        uniques[it.row()] = mtx.coeffRef(it.row(), cbidx) = it.value();
    }
    for (const std::pair<const std::vector<unsigned long long>, std::size_t>& umiset : mmaps_flat) {
        double factor;
        unsigned long long denom = 0;
        for (const unsigned long long &ftidx : umiset.first) {
            denom += uniques[ftidx];
        }
        if (denom == 0) {
            // Same as Uniform
            factor = 1.0;
            factor *= umiset.second;
            factor /= umiset.first.size();
            for (const unsigned long long &ftidx : umiset.first) {
                mtx.coeffRef(ftidx, cbidx) += factor;
            }
        } else {
            factor = 1.0 / denom;
            for (const unsigned long long &ftidx : umiset.first) {
                mtx.coeffRef(ftidx, cbidx) += uniques[ftidx] * umiset.second * factor;
            }
        }
    }
}

void DeduplicateWorker::handle_multimaps_EM(const int& cbidx, multimaps_mtx& mtx) {
    std::vector<double> uniques(numft);
    for (unique_mtx::InnerIterator it (unique_counts, cbidx); it; ++it) {
        uniques[it.row()] = it.value();
    }

    // Start with sum of unique and uniform
    std::vector<double> em_old = uniques, em(numft);
    for (const std::pair<const std::vector<unsigned long long>, std::size_t>& umiset : mmaps_flat) {
        double factor = 1.0;
        factor *= umiset.second;
        factor /= umiset.first.size();
        for (const unsigned long long &ftidx : umiset.first) {
            em_old[ftidx] += factor;
        }
    }
    // Use pointers to the two underlying vectors and interchange between them
    std::vector<double>* em1 = &em;
    std::vector<double>* em2 = &em_old;
    for (int iter = 0; iter < em_iter_max; ++iter) {
        std::vector<double>& vem1 = *em1, vem2 = *em2;

        // Clamp small values to 0
        std::replace_if(vem2.begin(), vem2.end(), [&](const double& x) {
            return x < em_small_thresh;
        }, 0.0);
        // Initialize with number of unique mappers
        vem1 = uniques;
        // For each pair (<ij>, N_i) do
        for (const std::pair<const std::vector<unsigned long long>, std::size_t>& umiset : mmaps_flat) {
            // Add N_i * normed(<C_ij>)
            double prev_sum = std::accumulate(umiset.first.cbegin(), umiset.first.cend(), 0.0, [&](const double& a, const unsigned long long& ftidx) -> double {
                return a + vem2[ftidx];
            });
            prev_sum = 1.0 / prev_sum;
            std::for_each(umiset.first.cbegin(), umiset.first.cend(), [&](const unsigned long long& ftidx) {
                vem1[ftidx] += vem2[ftidx] * umiset.second * prev_sum;
            });
        }
        // Test for convergence
        double maxval = std::inner_product(vem1.cbegin(), vem1.cend(), vem2.cbegin(), 0.0, [](const double& a, const double& x) -> double {
            return std::max(a, x);
        }, [](const double& a, const double &b) -> double {
            return std::abs(a - b);
        });
        std::swap(em1, em2); // Swap here to guarantee em2 points to last vector used
        if (maxval < em_max_diff) {
            break;
        }
    }
    // Populate nonzero entries from resulting vector
    std::vector<double>& vemlast = *em2;
    for (std::size_t ftidx = 0; ftidx < numft; ++ftidx) {
        if (vemlast[ftidx] != 0) {
            mtx.coeffRef(ftidx, cbidx) = vemlast[ftidx];
        }
    }
}

void DeduplicateWorker::handle_multimaps_Rescue(const int& cbidx, multimaps_mtx& mtx) {
    // Equivalent to one iteration of EM
    std::vector<double> uniques(numft);
    for (unique_mtx::InnerIterator it (unique_counts, cbidx); it; ++it) {
        uniques[it.row()] = mtx.coeffRef(it.row(), cbidx) = it.value();
    }
    for (const std::pair<const std::vector<unsigned long long>, std::size_t>& umiset : mmaps_flat) {
        double factor, uniffact = 1.0 / umiset.first.size();
        unsigned long long denom = 0;
        for (const unsigned long long &ftidx : umiset.first) {
            denom += uniques[ftidx];
        }
        denom += umiset.second;
        factor = 1.0 / denom;
        for (const unsigned long long &ftidx : umiset.first) {
            mtx.coeffRef(ftidx, cbidx) += umiset.second * (uniques[ftidx] + umiset.second * uniffact) * factor;
        }
    }
}

void DeduplicateWorker::handle_multimapped_reads(const int& cbidx) {
    if (multimapped_reads.empty()) {
        return;
    }
    // Convert the map<set<string>, size_t> to a vector<pair<vector<unsigned long long>, size_t>> flat array
    for (const std::pair<const std::set<std::string>, std::size_t>& pair : multimapped_reads) {
        std::vector<unsigned long long> &gene_idxs = mmaps_flat.emplace_back(0, pair.second).first;
        for (const std::string &name : pair.first) {
            gene_idxs.push_back(features[name].first);
        }
    }
    // Employ all requested multimap strategies
    for (auto& [strategy, matrix] : matrices) {
        switch (strategy) {
        case UNIQUE:
            std::unreachable();
            break;
        case UNIFORM:
            handle_multimaps_Uniform(cbidx, matrix);
            break;
        case PROPUNIQUE:
            handle_multimaps_PropUnique(cbidx, matrix);
            break;
        case EM:
            handle_multimaps_EM(cbidx, matrix);
            break;
        case RESCUE:
            handle_multimaps_Rescue(cbidx, matrix);
            break;
        }
    }
    // Reset for the next barcode
    mmaps_flat.clear();
    multimapped_reads.clear();
}

std::condition_variable time_printer_cv;
std::mutex time_printer_mtx;
unsigned long long time_printer_n_reads = 0;

void time_printer() {
    std::unique_lock lck{time_printer_mtx, std::defer_lock};
    // Update stderr every 30 seconds until the parent thread notifies
    while (time_printer_cv.wait_for(lck, 30s) == std::cv_status::timeout) {
        std::cerr << drugseq::put_time() << ": Processed " << time_printer_n_reads << " reads" << std::endl;
    }
}

int DeduplicateWorker::run() {
    constexpr int BCINFO_NO_INIT = -2;
    constexpr int BCINFO_NO_BARCODE = -1;
    constexpr int BCINFO_FIRST_VALID = 0;

    FilterKey umi;
    DedupStat result;
    BamTools::BamAlignment primary;
    std::string last_bc = "";
    std::pair<int, std::string> bcinfo = {BCINFO_NO_INIT, ""};
    {
    std::thread time_printer_task(time_printer);
    time_printer_task.detach();
    bool is_handling_multimappers = !matrices.empty();
    while (result = reader.get_read(umi, primary), result != DedupStat::FAIL_EOF && result != DedupStat::FAIL_HTS_ERROR) {
        if (last_bc != umi.cb) {
            if (!last_bc.empty()) {
                if (is_handling_multimappers) {
                    std::cerr << drugseq::put_time() << ": Handling multimaps for barcode " << last_bc << std::endl;
                    handle_multimapped_reads(bcinfo.first);
                }
                {
                    std::lock_guard guard{time_printer_mtx};
                    std::cerr << drugseq::put_time() << ": Finished barcode " << last_bc << " with " << time_printer_n_reads << " reads" << std::endl;
                }
            }
            last_bc = umi.cb;
            bcinfo = barcodes.at(last_bc);
            std::cerr << drugseq::put_time() << ": Begin processing UMIs for barcode " << last_bc << std::endl;
            {
                std::lock_guard guard{time_printer_mtx};
                time_printer_n_reads = 0;
            }
        }
        {
            std::lock_guard guard{time_printer_mtx};
            ++time_printer_n_reads;
        }
        std::vector<unsigned long long> &statrow = stats[bcinfo.first >= BCINFO_FIRST_VALID ? bcinfo.first : numbc];
        ++statrow[DedupStat::TOTAL_IN];
        ++statrow[result]; // QC review
        if (result == DedupStat::PASSED_QC) {
            if (umi.mx.size() > 1) {
                ++statrow[DedupStat::N_MULTIMAPPED];
            }
            if (bcinfo.first >= BCINFO_FIRST_VALID) {
                bool is_unique = umi.genes.size() == 1;
                if (is_unique) {
                    ++statrow[DedupStat::UNIQUE_COUNT];
                }
                std::optional<std::string> orig = detector.check_umi(umi);
                if (!orig) {
                    // Not a duplicate;
                    ++statrow[DedupStat::UMI_COUNT];
                    if (is_unique) {
                        ++statrow[DedupStat::UNIQUE_UMIS];
                        ++unique_counts.coeffRef(features[*umi.genes.cbegin()].first, bcinfo.first);
                    } else if (is_handling_multimappers) {
                        ++multimapped_reads[umi.genes];
                    }
                } else if (keep_bam) {
                    // Mark as duplicate and set tag DN = (qname of original read)
                    primary.SetIsDuplicate(true);
                    if (!primary.EditTag("DN", "Z", orig.value())) {
                        result = FAIL_HTS_ERROR;
                        break;
                    }
                }
            }
        } else if (keep_bam) {
            primary.SetIsFailedQC(true);
        }
        if (keep_bam) {
            if (!writer.SaveAlignment(primary)) {
                result = FAIL_HTS_ERROR;
                break;
            }
        }
    }
    // End the child thread by notifying the condition variable
    time_printer_cv.notify_all();
    }
    int ret = 0;
    if (result == FAIL_HTS_ERROR) {
        std::cerr << drugseq::put_time() << ": Fatal error with barcode " << last_bc << " after " << time_printer_n_reads << " reads" << std::endl;
        return -1;
    }
    std::cerr << drugseq::put_time() << ": Handling multimaps for barcode " << last_bc << std::endl;
    handle_multimapped_reads(bcinfo.first);
    // Since the child thread has already exited, we don't need to acquire the mutex here.
    std::cerr << drugseq::put_time() << ": Finished barcode " << last_bc << " with " << time_printer_n_reads << " reads" << std::endl;
    return 0;
}

void DeduplicateWorker::save_output(fs::path const& outdir) const {
    fs::create_directories(outdir);
    dump_barcodes(outdir / "barcodes.tsv");
    dump_features(outdir / "features.tsv");
    dump_stats(outdir / "dedup_stats.tsv");
    DeduplicateWorker::dump_matrix(outdir / "matrix.mtx", unique_counts);
    for (const std::pair<MultimapStrategy, multimaps_mtx>& pair : matrices) {
        fs::path mtx_fname = outdir;
        bool valid = true;
        switch (pair.first) {
        case MultimapStrategy::UNIFORM:
            mtx_fname /= "UniqueAndMult-Uniform.mtx";
            break;
        case MultimapStrategy::PROPUNIQUE:
            mtx_fname /= "UniqueAndMult-PropUnique.mtx";
            break;
        case MultimapStrategy::EM:
            mtx_fname /= "UniqueAndMult-EM.mtx";
            break;
        case MultimapStrategy::RESCUE:
            mtx_fname /= "UniqueAndMult-Rescue.mtx";
            break;
        default:
            valid = false;
            break;
        }
        if (valid) {
            DeduplicateWorker::dump_matrix(mtx_fname, pair.second);
        }
    }
}

void DeduplicateWorker::dump_barcodes(const fs::path& fname) const {
    std::ofstream outfile(fname);
    std::vector<std::string> bcvec(numbc);
    for (const std::pair<const std::string, std::pair<int, std::string>> bcinfo : barcodes) {
        if (bcinfo.second.first >= 0) {
            bcvec[bcinfo.second.first] = bcinfo.first;
        }
    }
    for (const std::string &bc : bcvec) {
        outfile << bc << "\n";
    }
}

void DeduplicateWorker::dump_features(const fs::path& fname) const {
    std::ofstream outfile(fname);
    std::vector<std::pair<std::string, std::string>> ftvec(numft);
    for (const std::pair<const std::string, std::pair<int, std::string>> ftinfo : features) {
        ftvec[ftinfo.second.first] = {ftinfo.first, ftinfo.second.second};
    }
    for (const std::pair<std::string, std::string> &ft : ftvec) {
        outfile << ft.first << "\t" << ft.second << "\t" << "Gene Expression\n";
    }
}

void DeduplicateWorker::dump_stats(const fs::path& fname) const {
    std::ofstream outfile(fname);
    std::vector<std::string> bcvec(numbc + 1);
    for (const std::pair<const std::string, std::pair<int, std::string>> bcinfo : barcodes) {
        bcvec[bcinfo.second.first >= 0 ? bcinfo.second.first : numbc] = bcinfo.second.second;
    }
    outfile << "\tTOTAL\tDISCARD_QCFAIL\tDISCARD_UNMAPPED\tDISCARD_NOFEATURE\tPASSED_QC\tN_MULTIMAPPED\tUMI_COUNTS\tUNIQUE_COUNTS\tUNIQUE_UMIS\n";
    for (std::size_t i = 0; i <= numbc; ++i) {
        outfile << bcvec[i];
        for (const unsigned long long &s : stats[i]) {
            outfile << "\t" << s;
        }
        outfile << "\n";
    }
}
