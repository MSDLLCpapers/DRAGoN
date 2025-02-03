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
#include <iostream>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <drugseq.hpp>
#include "DemultiplexWorker.hpp"

constexpr std::string_view usage =
"usage: demultiplex [-h] [-v] [--bclen BCLEN] [--umistart UMISTART] [--umilen UMILEN] [--split SPLIT] [--mismatch MISMATCH] [--min-length MIN_LENGTH] [--min-qual MIN_QUAL] [--min-qual-bases MIN_QUAL_BASES] [--max-AT MAX_AT] [--max-GC MAX_GC]\n"
"                   [--min-homop-len MIN_HOMOP_LEN] [--max-homop-mmatch MAX_HOMOP_MMATCH] [--max-bc-ns MAX_BC_NS] [--check-bcs] [--cpus CPUS] [--adapter ADAPTER]\n"
"                   R1 R2 barcodes prefix\n";
constexpr std::string_view more_usage =
"\n"
"positional arguments:\n"
"  R1                    Path to read-1 fastq file (barcode+UMI)\n"
"  R2                    Path to read-2 fastq file (cDNA)\n"
"  barcodes              Path to TSV file defining well names and barcode sequences\n"
"  prefix                Output file prefix\n"
"\n"
"options:\n"
"  -h, --help            show this help message and exit\n"
"  -v, --version         show the program version and exit\n"
"  --bclen BCLEN         Number of bases in the well barcode (default: 16)\n"
"  --umistart UMISTART   Start position of the UMI in the R1 read (default: 17)\n"
"  --umilen UMILEN       Number of bases in the UMI sequence (default: 10)\n"
"  --split SPLIT         Maximum number of BAM files to write (default: one per well)\n"
"  --mismatch MISMATCH   Maximum number of mismatches (substitutions) between observed and defined barcode sequences (default: 1)\n"
"  --min-length MIN_LENGTH\n"
"                        Minimum cDNA length to pass QC, before or after trimming (default: 20)\n"
"  --min-qual MIN_QUAL   Phred quality score threshold for UMI filtering and cDNA trimming (default: 20)\n"
"  --min-qual-bases MIN_QUAL_BASES\n"
"                        Minimum number of UMI bases with phred quality score passing --min-qual (default: 6)\n"
"  --max-AT MAX_AT       Maximum fraction of cDNA bases allowed to be A or T per read (default: 0.9)\n"
"  --max-GC MAX_GC       Maximum fraction of cDNA bases allowed to be G or C per read (default: 0.9)\n"
"  --min-homop-len MIN_HOMOP_LEN\n"
"                        Trim the cDNA when I encounter a run of at least this many of the same nucleotide (homopolymer) (default: 10)\n"
"  --max-homop-mmatch MAX_HOMOP_MMATCH\n"
"                        Allow at most this many substitutions when detecting homopolymers (default: 0)\n"
"  --max-bc-ns MAX_BC_NS\n"
"                        Maximum number of N base calls allowed in the barcode sequence (default: 1)\n"
"  --check-bcs           If set, will make sure that the barcodes are unique, and the list of barcodes is compatible with the --mismatch setting.\n"
"  --cpus CPUS           Number of BAM IO threads\n"
"  --adapter ADAPTER     Sequencing adapter to trim from the 3' end of R2\n"
"  --discard-failed      If set, will excluded QC-failed reads from the BAM files. Otherwise, will output a separate BAM file with QC-fail reads.";

typedef int opt_parse_t;
constexpr opt_parse_t opt_no = 0, opt_yes = 1, opt_inval = 2;

opt_parse_t get_arg(std::vector<std::string>::const_iterator& it, const std::vector<std::string>::const_iterator& end, std::string const& optstr, std::string& out) {
    if (*it != optstr) {
        return opt_no;
    }
    if (++it == end || it->front() == '-') {
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

opt_parse_t get_arg(std::vector<std::string>::const_iterator& it, const std::vector<std::string>::const_iterator& end, std::string const& optstr, int& out) {
    if (*it != optstr) {
        return opt_no;
    }
    if (++it == end || it->front() == '-') {
        std::cerr << usage << "\nERR: missing value for " << optstr << std::endl;
        return opt_inval;
    }
    try {
        out = std::stoi(*it);
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
    if (++it == end || it->front() == '-') {
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

int main(int argc, char ** argv) {
    int umilen = 10;
    int split = 0;
    int max_mismatch = 1;
    int max_n = 1;
    int min_len = 20;
    int min_qual = 20;
    int phred = 33;
    int umi_qual_bases = 6;
    int homop_len = 10;
    double max_AT = 0.9;
    double max_GC = 0.9;
    int bamthreads = 1;
    int dummy;
    bool write_qcfail = true;
    bool debug = false;
    std::string adapter = "";
    opt_parse_t opt_parse_result;
    std::vector<std::string> args(argv, argv + argc);
    std::vector<std::string> posargs {}; // R1 R2 barcodes prefix

    for (auto it = std::next(args.cbegin()); it != args.cend(); ++it) {
        const std::string& opt = *it;
        if (opt == "-h" || opt == "--help") {
            std::cout << usage << more_usage << "\n";
            return 0;
        } else if (opt == "-v" || opt == "--version") {
            std::cout << "deduplicate (DRAGoN v" << DRUGSEQ_VERSION_STR << ")" << std::endl;
            return 0;
        } else if (opt == "--debug") {
            debug = true;
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--bclen", dummy)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--umistart", dummy)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--umilen", umilen)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--split", split)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--mismatch", max_mismatch)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--min-length", min_len)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--min-qual", min_qual)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--min-qual-bases", umi_qual_bases)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--max-AT", max_AT)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--max-GC", max_GC)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--min-homop-len", homop_len)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--max-homop-mmatch", dummy)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--max-bc-ns", max_n)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if (opt == "--check-bcs") {
            // legacy
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--cpus", bamthreads)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if ((opt_parse_result = get_arg(it, args.cend(), "--adapter", adapter)) != opt_no) {
            if (opt_parse_result == opt_inval) {
                return 1;
            }
        } else if (opt == "--discard-failed") {
            write_qcfail = false;
        } else if (opt.front() == '-') {
            std::cerr << usage << "\nERR: Unrecognized option: " << opt << "\n";
            return 1;
        } else {
            posargs.push_back(opt);
        }
    }
    if (posargs.size() < 4) {
        const std::string posarg_name[4] {
            "R1",
            "R2",
            "barcodes",
            "prefix"
        };
        std::cerr << usage << "\nERR: Missing required argument: " << posarg_name[posargs.size()] << "\n";
        return 1;
    } else if (posargs.size() > 4) {
        std::cerr << usage << "\nERR: Unrecognized positional argument (first one: " << posargs[4] << ")\n";
        return 1;
    }

    // Check for invalid argument values
    if (umilen < 1) {
        std::cerr << usage << "\nERR: Minimum allowed value for --umilen is 1\n";
        return 1;
    }
    if (max_mismatch < 0) {
        std::cerr << usage << "\nERR: Minimum allowed value for --mismatch is 0\n";
        return 1;
    }
    if (min_len < 0) {
        std::cerr << usage << "\nERR: Minimum allowed value for --min-length is 0\n";
        return 1;
    }
    if (min_len < 20) {
        std::cerr << "WARN: Setting --min-length less than 20 is not advised\n";
    }
    if (min_qual < 0) {
        std::cerr << usage << "\nERR: Minimum allowed value for --min-qual is 0\n";
        return 1;
    }
    if (umi_qual_bases < 0 || umi_qual_bases > umilen) {
        std::cerr << usage << "\nERR: Allowed values for --min-qual-bases is between 0 and " << umilen << "\n";
        return 1;
    }
    if (homop_len < 0) {
        std::cerr << usage << "\nERR: Minimum allowed value for --min-homop-len is 0\n";
        return 1;
    }
    if (max_AT < 0 || max_AT > 1) {
        std::cerr << usage << "\nERR: Value for --max-AT must be between 0 and 1\n";
        return 1;
    }
    if (max_GC < 0 || max_GC > 1) {
        std::cerr << usage << "\nERR: Value for --max-GC must be between 0 and 1\n";
        return 1;
    }

    boost::log::trivial::severity_level level = debug ? boost::log::trivial::debug : boost::log::trivial::info;
    boost::log::core::get()->set_filter (
        boost::log::trivial::severity >= level
    );

    // Do the thing
    std::string bamfpref = posargs.at(3);
    std::ofstream json (posargs.at(3) + ".demux.json");
    try {
        MainWorker worker {
            posargs.at(2),
            bamfpref,
            args,
            umilen,
            split,
            max_mismatch,
            max_n,
            min_len,
            min_qual,
            phred,
            umi_qual_bases,
            homop_len,
            max_AT,
            max_GC,
            bamthreads,
            adapter,
            write_qcfail
        };
        if (worker.run(posargs.at(0), posargs.at(1)) < 0) {
            std::cerr << "ERR: Demultiplexing failed" << std::endl;
            return 1;
        }
        worker.dump_stats(json);
    } catch (const std::exception& e) {
        std::cerr << "ERR: Demultiplexing failed: " << e.what() << std::endl;
        return 1;
    }
}
