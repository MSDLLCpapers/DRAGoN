/*
 * Written by Scott Norton
 * 15-Dec-2022
 */

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <unordered_set>
#include "MultimapHandlers.hpp"
#include "DeduplicateWorker.hpp"

constexpr std::string_view usage =
"usage: deduplicate [-h] [-v] [--cpus CPUS] [--umilen UMILEN]\n"
"                   [--max-distance MAX_DISTANCE]\n"
"                   [--max-lookback MAX_LOOKBACK]\n"
"                   [--mmap-strategies {Unique,Uniform,PropUnique,EM,Rescue} [{Unique,Uniform,PropUnique,EM,Rescue} ...]]\n"
"                   [--em-iter-max EM_ITER_MAX]\n"
"                   [--em-max-diff EM_MAX_DIFF]\n"
"                   [--em-small-thresh EM_SMALL_THRESH] [--debug] [--keep-bam]\n"
"                   bamfile barcodes annofile outdir\n";
constexpr std::string_view more_usage =
"\n"
"positional arguments:\n"
"  bamfile               Annotated BAM file\n"
"  barcodes              TSV file mapping barcode sequences to well names and\n"
"                        coordinates\n"
"  annofile              GTF file defining quantified features\n"
"  outdir                Output directory\n"
"\n"
"options:\n"
"  -h, --help            show this help message and exit\n"
"  -v, --version         show the program version and exit\n"
"  --cpus CPUS           Number of threads for BAM IO (default: 1)\n"
"  --umilen UMILEN       UMI sequence length (default: 10)\n"
"  --max-distance MAX_DISTANCE\n"
"                        Maximum edit distance between two UMIs to consider\n"
"                        them as duplicates (default: 2)\n"
"  --max-lookback MAX_LOOKBACK\n"
"                        Maximum genomic distance between two UMIs to test them\n"
"                        as duplicates (default: 10)\n"
"  --mmap-strategies {Unique,Uniform,PropUnique,EM,Rescue} [{Unique,Uniform,PropUnique,EM,Rescue} ...]\n"
"                        Strategy(ies) for handling reads assigned to multiple\n"
"                        features (default: ['Unique'])\n"
"  --em-iter-max EM_ITER_MAX\n"
"                        If --ambiguous contains EM, the maximum number of\n"
"                        iterations of EM (default: 100)\n"
"  --em-max-diff EM_MAX_DIFF\n"
"                        If --ambiguous contains EM, the tolerance for\n"
"                        convergence, measured by the largest absolute change\n"
"                        in counts for any gene in the current well (default:\n"
"                        0.01)\n"
"  --em-small-thresh EM_SMALL_THRESH\n"
"                        If --ambiguous contains EM, clamp counts smaller than\n"
"                        this to 0 (default: 0.01)\n"
"  --debug               Increase logging verbosity\n";

typedef int opt_parse_t;
constexpr opt_parse_t opt_no = 0, opt_yes = 1, opt_inval = 2;

opt_parse_t get_arg(std::vector<std::string>::const_iterator& it, const std::vector<std::string>::const_iterator& end, std::string const& optstr, std::string& out) {
    if (*it != optstr) {
        return opt_no;
    }
    if (++it == end || (*it)[0] == '-') {
        std::cerr << usage << "\nERR: missing value for " << optstr << std::endl;
        return opt_inval;
    }
    try {
        out = *it;
    } catch (const std::exception& e) {
        std::cerr << usage << "\nERR: invalid value for " << optstr << ": " << e.what() << std::endl;
        return opt_inval;
    }
    return opt_yes;
}

opt_parse_t get_arg(std::vector<std::string>::const_iterator& it, const std::vector<std::string>::const_iterator& end, std::string const& optstr, unsigned long long& out) {
    if (*it != optstr) {
        return opt_no;
    }
    if (++it == end || (*it)[0] == '-') {
        std::cerr << usage << "\nERR: missing value for " << optstr << std::endl;
        return opt_inval;
    }
    try {
        out = std::stoull(*it);
    } catch (const std::exception& e) {
        std::cerr << usage << "\nERR: invalid value for " << optstr << ": " << e.what() << std::endl;
        return opt_inval;
    }
    return opt_yes;
}

opt_parse_t get_arg(std::vector<std::string>::const_iterator& it, const std::vector<std::string>::const_iterator& end, std::string const& optstr, double& out) {
    if (*it != optstr) {
        return opt_no;
    }
    if (++it == end || (*it)[0] == '-') {
        std::cerr << usage << "\nERR: missing value for " << optstr << std::endl;
        return opt_inval;
    }
    try {
        out = std::stod(*it);
    } catch (const std::exception& e) {
        std::cerr << usage << "\nERR: invalid value for " << optstr << ": " << e.what() << std::endl;
        return opt_inval;
    }
    return opt_yes;
}

opt_parse_t get_arg(std::vector<std::string>::const_iterator& it, const std::vector<std::string>::const_iterator& end, std::string const& optstr, std::unordered_set<MultimapStrategy>& out) {
    if (*it != optstr) {
        return opt_no;
    }
    if (++it == end || (*it)[0] == '-') {
        std::cerr << usage << "\nERR: missing value for " << optstr << std::endl;
        return opt_inval;
    }
    try {
        do {
            std::string const& x = *it++;
            if (x == "Unique") {
                // UNIQUE is implied
            } else if (x == "Uniform") {
                out.insert(MultimapStrategy::UNIFORM);
            } else if (x == "PropUnique") {
                out.insert(MultimapStrategy::PROPUNIQUE);
            } else if (x == "EM") {
                out.insert(MultimapStrategy::EM);
            } else if (x == "Rescue") {
                out.insert(MultimapStrategy::RESCUE);
            } else {
                std::cerr << usage << "\nERR: Invalid MultimapStrategy value: \"" << x << "\"\n  Choose from Unique, Uniform, PropUnique, EM, Rescue\n";
                return 1;
            }
        } while (it != end && (*it)[0] != '-');
    } catch (const std::exception& e) {
        std::cerr << usage << "\nERR: invalid value for " << optstr << ": " << e.what() << std::endl;
        return opt_inval;
    }
    --it;
    return opt_yes;
}

int main(int argc, char ** argv) {
    std::cerr << drugseq::put_time() << ": Begin deduplicate workflow" << std::endl;
    unsigned long long threads = 1;
    std::unordered_set<MultimapStrategy> mmap_strategies;
    unsigned long long umilen = 12;
    unsigned long long max_distance = 2;
    unsigned long long max_lookback = 10;
    unsigned long long em_iter_max = 100;
    double em_max_diff = 0.01;
    double em_small_thresh = 0.01;
    bool keep_bam = false;
    opt_parse_t opt_parse_result;

    std::vector<std::string> cli(argv, argv + argc);
    std::vector<std::string> posargs; // bamfile barcode annofile outdir

    for (auto it = std::next(cli.cbegin()); it != cli.cend(); ++it) {
        const std::string& opt = *it;
        if (opt == "-h" || opt == "--help") {
            std::cout << usage << more_usage << std::endl;
            return 0;
        } else if (opt == "-v" || opt == "--version") {
            std::cout << "deduplicate (DRAGoN v" << DRUGSEQ_VERSION_STR << ")" << std::endl;
            return 0;
        } else if ((opt_parse_result = get_arg(it, cli.cend(), "--cpus", threads)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, cli.cend(), "--umilen", umilen)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, cli.cend(), "--max-distance", max_distance)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, cli.cend(), "--max-lookback", max_lookback)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, cli.cend(), "--mmap-strategies", mmap_strategies)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, cli.cend(), "--em-iter-max", em_iter_max)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, cli.cend(), "--em-max-diff", em_max_diff)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, cli.cend(), "--em-small-thresh", em_small_thresh)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if (opt == "--keep-bam") {
            keep_bam = true;
        } else if (opt[0] == '-') {
            std::cerr << usage << "\nERR: Unrecognized option flag: \"" << opt << "\"\n";
            return 1;
        } else {
            posargs.push_back(opt);
        }
    }

    if (posargs.size() < 4) {
        const std::string argnames[] {"bamfile", "barcodes", "annofile", "outdir"};
        std::cerr << usage << "\nERR: Missing required positional opt: " << argnames[posargs.size()] << "\n";
        return 1;
    }
    if (posargs.size() > 4) {
        std::cerr << usage << "\nERR: Extra unrecognized positional opt (first one: " << posargs[4] << ")\n";
        return 1;
    }

    // Check for invalid argument values
    if (umilen < 1) {
        std::cerr << usage << "\nERR: Minimum allowed value for --umilen is 1\n";
        return 1;
    }
    if (max_distance > umilen / 2) {
        std::cerr << usage << "\nERR: When setting --umilen " << umilen << ", maximum allowed value for --max-distance is " << umilen / 2 << "\n";
        return 1;
    }

    try {
        DeduplicateWorker worker(
            cli,
            posargs[0],
            threads,
            posargs[1],
            posargs[2],
            mmap_strategies,
            umilen,
            max_distance,
            max_lookback,
            em_iter_max,
            em_max_diff,
            em_small_thresh,
            keep_bam
        );
        int result = worker.run();
        if (result != 0) {
            std::cerr << "ERR: Deduplication failed\n";
            return result;
        }
        worker.save_output(posargs[3]);
        std::cerr << drugseq::put_time() << ": Finished deduplicate workflow!" << std::endl;
    } catch (std::exception& e) {
        std::cerr << "ERR: Error deduplicating: " << e.what() << std::endl;
        return 1;
    }
}
