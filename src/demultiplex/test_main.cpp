// Written by Scott Norton
// 15-Feb-2023

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#include <boost/test/unit_test.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <kseq++/seqio.hpp>
#include <bamtools/api/BamReader.h>
#include <bamtools/api/BamAlignment.h>
#include "DemultiplexWorker.hpp"

template<class T>
T get_tag_value(BamTools::BamAlignment const& read, std::string const& tag) {
    T ret;
    BOOST_CHECK(read.GetTag<T>(tag, ret));
    return ret;
}

namespace fs = std::filesystem;

fs::path make_temp_filename(const std::string& suffix) {
    static char path[PATH_MAX+1] {};
    int fd = mkstemps(const_cast<char*>((fs::temp_directory_path() / ("tmpXXXXXX" + suffix)).c_str()), suffix.length());
    if (fd == -1) {
        throw std::runtime_error("io error");
    }
    std::string link = "/proc/self/fd/" + std::to_string(fd);
    if (readlink(link.c_str(), path, PATH_MAX) == -1) {
        throw std::runtime_error("io error");
    }
    return path;
}

struct run_demultiplex_single {
    size_t cdna_n_good_bases = -1ull;
    size_t umi_n_good_bases = -1ull;
    std::string adapter = "";
    bool expected_failure = false;

    BamTools::BamAlignment& operator()(
        BamTools::BamAlignment& read,
        std::string const& cdnaseq,
        std::string const& bcseq,
        std::string const& umiseq
    ) const {
        boost::log::core::get()->set_filter (
            boost::log::trivial::severity >= boost::log::trivial::debug
        );

        fs::path barcodes_fn = make_temp_filename(".txt");
        fs::path bam_stem = make_temp_filename("");
        fs::path fq1_fn = make_temp_filename(".fastq.gz");
        fs::path fq2_fn = make_temp_filename(".fastq.gz");
        {
            std::ofstream barcodes_fp(barcodes_fn);
            barcodes_fp <<
            "AAACGGCGAGCCCTTT\tPool1_A1__Media1\t1\t1\n"
            "AGACTTTCAATCACTA\tPool1_B1__293T\t1\t2\n"
            "CCTACCATCAAGTCCA\tPool1_C1__293T\t1\t3\n"
            "CCGAGGCACCGAGAAC\tPool1_D1__293T\t1\t4\n"
            "GGCGTTTCGTGAGAAG\tPool1_E1__293T_DMSO\t1\t5\n"
            "GATGCAGTGTAGTGCG\tPool1_F1__293T_DMSO\t1\t6\n"
            "TTGTAAGATCTTAGGC\tPool1_G1__293T_DMSO\t1\t7\n"
            "TTCTACAGTGCTCTGT\tPool1_H1__DMSO\t1\t8\n"
            "CTGCAGCCCATGGCAG\tPool1_A2__Media1\t2\t1\n"
            "CGTAGGTTCTTAGTCA\tPool1_B2__293T\t2\t2\n";
        }
        {
            klibpp::SeqStreamOut r1fp(fq1_fn.c_str());
            std::string qual_bases (bcseq.size() + std::min(umi_n_good_bases, umiseq.size()), 'F');
            if (umi_n_good_bases < umiseq.size()) {
                qual_bases += std::string(umiseq.size() - umi_n_good_bases, '*');
            }
            r1fp << (klibpp::KSeq {"TestRead", "R1", bcseq + umiseq, qual_bases});
        }
        {
            klibpp::SeqStreamOut r2fp(fq2_fn.c_str());
            std::string qual_bases (std::min(cdna_n_good_bases, cdnaseq.size()), 'F');
            if (cdna_n_good_bases < cdnaseq.size()) {
                qual_bases += std::string(cdnaseq.size() - cdna_n_good_bases, '*');
                BOOST_TEST_MESSAGE("CDNA qual bases = " << qual_bases);
            }
            r2fp << (klibpp::KSeq {"TestRead", "R2", cdnaseq, qual_bases});
        }

        std::vector<std::string> cli {"demultiplex", fq1_fn, fq2_fn, barcodes_fn, bam_stem, "--umilen", std::to_string(umiseq.size())};
        {
            MainWorker worker {
                barcodes_fn,
                bam_stem.string(),
                cli,
                static_cast<int>(umiseq.size()),
                0,
                1,
                1,
                20,
                20,
                33,
                6,
                10,
                0.9,
                0.9,
                1,
                adapter,
                true
            };
            BOOST_REQUIRE_EQUAL((worker.run(fq1_fn, fq2_fn)), 0);
            {
                std::ofstream stats_file (bam_stem.string() + ".json");
                worker.dump_stats(stats_file);
            }
        }
        BOOST_TEST_PASSPOINT();

        {
            fs::path path;
            path = bam_stem.string() + (expected_failure ? "_qcFail.bam" : ".bam");
            BOOST_TEST_REQUIRE(fs::exists(path));
            BamTools::BamReader samfile {};
            BOOST_TEST_REQUIRE(!samfile.IsOpen());
            BOOST_TEST_REQUIRE(samfile.Open(path.string()));
            BOOST_TEST_REQUIRE(samfile.GetNextAlignment(read));
            samfile.Close();
        }
        return read;
    }
};

BOOST_AUTO_TEST_CASE(test_demultiplex) {
    BamTools::BamAlignment read {};
    run_demultiplex_single{}(
        read,
        "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA",
        "TTGTAAGATCTTAGGC",
        "CATACGACTGCA"
    );
    BOOST_CHECK_EQUAL(read.Name, "TestRead");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CR"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UR"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CY"), "FFFFFFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UY"), "FFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CB"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UB"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "BI"), "0006");
    BOOST_CHECK_EQUAL(get_tag_value<int>(read, "BM"), 0);
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "QC"), "PASS");
    BOOST_CHECK_EQUAL(read.QueryBases, "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA");
    BOOST_CHECK_EQUAL(read.Qualities, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
    BOOST_CHECK(!read.IsMapped());
}

BOOST_AUTO_TEST_CASE(test_demultiplex_1mm) {
    BamTools::BamAlignment read {};
    run_demultiplex_single{}(
        read,
        "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA",
        "TTGTCAGATCTTAGGC",
        "CATACGACTGCA"
    );
    BOOST_CHECK_EQUAL(read.Name, "TestRead");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CR"), "TTGTCAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UR"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CY"), "FFFFFFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UY"), "FFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CB"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UB"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "BI"), "0006");
    BOOST_CHECK_EQUAL(get_tag_value<int>(read, "BM"), 1);
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "QC"), "PASS");
    BOOST_CHECK_EQUAL(read.QueryBases, "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA");
    BOOST_CHECK_EQUAL(read.Qualities, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
    BOOST_CHECK(!read.IsMapped());
}

BOOST_AUTO_TEST_CASE(test_demultiplex_noBC) {
    BamTools::BamAlignment read {};
    run_demultiplex_single{}(
        read,
        "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA",
        "TTGTCAGATCTTAGGG",
        "CATACGACTGCA"
    );
    BOOST_CHECK_EQUAL(read.Name, "TestRead");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CR"), "TTGTCAGATCTTAGGG");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UR"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CY"), "FFFFFFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UY"), "FFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CB"), "NNNNNNNNNNNNNNNN");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UB"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "BI"), "noBC");
    BOOST_CHECK_EQUAL(get_tag_value<int>(read, "BM"), -1);
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "QC"), "PASS");
    BOOST_CHECK_EQUAL(read.QueryBases, "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA");
    BOOST_CHECK_EQUAL(read.Qualities, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
    BOOST_CHECK(!read.IsMapped());
}

BOOST_AUTO_TEST_CASE(test_demultiplex_fail_homopolymer) {
    BamTools::BamAlignment read {};
    run_demultiplex_single{
        .expected_failure = true
    }(
        read,
        "AGCTACGTAGTAGCAGCAAAAAAAAAAAATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAATAGCATAGCTAG",
        "TTGTAAGATCTTAGGC",
        "CATACGACTGCA"
    );
    BOOST_CHECK_EQUAL(read.Name, "TestRead");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CR"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UR"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CY"), "FFFFFFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UY"), "FFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CB"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UB"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "BI"), "0006");
    BOOST_CHECK_EQUAL(get_tag_value<int>(read, "BM"), 0);
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "QC"), "FAIL_HOMOPOLYMER");
    BOOST_CHECK_EQUAL(read.QueryBases, "AGCTACGTAGTAGCAGCAAAAAAAAAAAATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAATAGCATAGCTAG");
    BOOST_CHECK_EQUAL(read.Qualities, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
    BOOST_CHECK(!read.IsMapped());
    BOOST_CHECK(read.IsFailedQC());
}

BOOST_AUTO_TEST_CASE(test_demultiplex_pass_homopolymer) {
    BamTools::BamAlignment read {};
    run_demultiplex_single{}(
        read,
        "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTGGGGGGGGGGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA",
        "TTGTAAGATCTTAGGC",
        "CATACGACTGCA"
    );
    BOOST_CHECK_EQUAL(read.Name, "TestRead");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CR"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UR"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CY"), "FFFFFFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UY"), "FFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CB"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UB"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "BI"), "0006");
    BOOST_CHECK_EQUAL(get_tag_value<int>(read, "BM"), 0);
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "QC"), "PASS");
    BOOST_CHECK_EQUAL(read.QueryBases, "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGT");
    BOOST_CHECK_EQUAL(read.Qualities, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
    BOOST_CHECK(!read.IsMapped());
}

BOOST_AUTO_TEST_CASE(test_demultiplex_fail_length) {
    BamTools::BamAlignment read {};
    run_demultiplex_single{
        .expected_failure=true
    }(
        read,
        "AGCTACGTAGTAGCAGCTA",
        "TTGTAAGATCTTAGGC",
        "CATACGACTGCA"
    );
    BOOST_CHECK_EQUAL(read.Name, "TestRead");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CR"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UR"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CY"), "FFFFFFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UY"), "FFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CB"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UB"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "BI"), "0006");
    BOOST_CHECK_EQUAL(get_tag_value<int>(read, "BM"), 0);
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "QC"), "FAIL_LENGTH");
    BOOST_CHECK_EQUAL(read.QueryBases, "AGCTACGTAGTAGCAGCTA");
    BOOST_CHECK_EQUAL(read.Qualities, "FFFFFFFFFFFFFFFFFFF");
    BOOST_CHECK(!read.IsMapped());
    BOOST_CHECK(read.IsFailedQC());
}

BOOST_AUTO_TEST_CASE(test_demultiplex_fail_trimming) {
    BamTools::BamAlignment read {};
    run_demultiplex_single{
        .cdna_n_good_bases = 10,
        .expected_failure = true
    }(
        read,
        "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA",
        "TTGTAAGATCTTAGGC",
        "CATACGACTGCA"
    );
    BOOST_CHECK_EQUAL(read.Name, "TestRead");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CR"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UR"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CY"), "FFFFFFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UY"), "FFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CB"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UB"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "BI"), "0006");
    BOOST_CHECK_EQUAL(get_tag_value<int>(read, "BM"), 0);
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "QC"), "FAIL_QUALITY_TRIMMING");
    BOOST_CHECK_EQUAL(read.QueryBases, "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA");
    BOOST_CHECK_EQUAL(read.Qualities, "FFFFFFFFFF********************************************************************************");
    BOOST_CHECK(!read.IsMapped());
    BOOST_CHECK(read.IsFailedQC());
}

BOOST_AUTO_TEST_CASE(test_demultiplex_pass_trimming) {
    BamTools::BamAlignment read {};
    run_demultiplex_single{
        .cdna_n_good_bases = 50
    }(
        read,
        "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA",
        "TTGTAAGATCTTAGGC",
        "CATACGACTGCA"
    );
    BOOST_CHECK_EQUAL(read.Name, "TestRead");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CR"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UR"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CY"), "FFFFFFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UY"), "FFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CB"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UB"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "BI"), "0006");
    BOOST_CHECK_EQUAL(get_tag_value<int>(read, "BM"), 0);
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "QC"), "PASS");
    BOOST_CHECK_EQUAL(read.QueryBases, "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTA");
    BOOST_CHECK_EQUAL(read.Qualities, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
    BOOST_CHECK(!read.IsMapped());
}

BOOST_AUTO_TEST_CASE(test_demultiplex_fail_umi_quality) {
    BamTools::BamAlignment read {};
    run_demultiplex_single{
        .umi_n_good_bases = 5,
        .expected_failure = true
    }(
        read,
        "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA",
        "TTGTAAGATCTTAGGC",
        "CATACGACTGCA"
    );
    BOOST_CHECK_EQUAL(read.Name, "TestRead");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CR"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UR"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CY"), "FFFFFFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UY"), "FFFFF*******");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CB"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UB"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "BI"), "0006");
    BOOST_CHECK_EQUAL(get_tag_value<int>(read, "BM"), 0);
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "QC"), "FAIL_UMI_QUALITY");
    BOOST_CHECK_EQUAL(read.QueryBases, "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA");
    BOOST_CHECK_EQUAL(read.Qualities, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
    BOOST_CHECK(!read.IsMapped());
    BOOST_CHECK(read.IsFailedQC());
}

BOOST_AUTO_TEST_CASE(test_demultiplex_fail_adapter_trimming) {
    BamTools::BamAlignment read {};
    run_demultiplex_single{
        .adapter = "CAGCTAGCATAGCTAG",
        .expected_failure = true
    }(
        read,
        "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA",
        "TTGTAAGATCTTAGGC",
        "CATACGACTGCA"
    );
    BOOST_CHECK_EQUAL(read.Name, "TestRead");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CR"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UR"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CY"), "FFFFFFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UY"), "FFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CB"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UB"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "BI"), "0006");
    BOOST_CHECK_EQUAL(get_tag_value<int>(read, "BM"), 0);
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "QC"), "FAIL_ADAPTER_TRIMMING");
    BOOST_CHECK_EQUAL(read.QueryBases, "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA");
    BOOST_CHECK_EQUAL(read.Qualities, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
    BOOST_CHECK(!read.IsMapped());
    BOOST_CHECK(read.IsFailedQC());
}

BOOST_AUTO_TEST_CASE(test_demultiplex_pass_adapter_trimming) {
    BamTools::BamAlignment read {};
    run_demultiplex_single{
        .adapter = "CTAGTCGATCGATCGT"
    }(
        read,
        "AGCTACGTAGTAGCAGCTAGCATAGCTAGCTAGTCGATCGATCGTACGTAGCTAGCATGCTAGCTAGCTAGCTAGCTAATCTACGTACAA",
        "TTGTAAGATCTTAGGC",
        "CATACGACTGCA"
    );
    BOOST_CHECK_EQUAL(read.Name, "TestRead");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CR"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UR"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CY"), "FFFFFFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UY"), "FFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CB"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UB"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "BI"), "0006");
    BOOST_CHECK_EQUAL(get_tag_value<int>(read, "BM"), 0);
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "QC"), "PASS");
    BOOST_CHECK_EQUAL(read.QueryBases, "AGCTACGTAGTAGCAGCTAGCATAGCTA");
    BOOST_CHECK_EQUAL(read.Qualities, "FFFFFFFFFFFFFFFFFFFFFFFFFFFF");
    BOOST_CHECK(!read.IsMapped());
}

BOOST_AUTO_TEST_CASE(test_demultiplex_fail_at_content) {
    BamTools::BamAlignment read {};
    run_demultiplex_single{
        .expected_failure = true
    }(
        read,
        "AGCTACGTATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT",
        "TTGTAAGATCTTAGGC",
        "CATACGACTGCA"
    );
    BOOST_CHECK_EQUAL(read.Name, "TestRead");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CR"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UR"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CY"), "FFFFFFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UY"), "FFFFFFFFFFFF");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "CB"), "TTGTAAGATCTTAGGC");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "UB"), "CATACGACTGCA");
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "BI"), "0006");
    BOOST_CHECK_EQUAL(get_tag_value<int>(read, "BM"), 0);
    BOOST_CHECK_EQUAL(get_tag_value<std::string>(read, "QC"), "FAIL_AT_CONTENT");
    BOOST_CHECK_EQUAL(read.QueryBases, "AGCTACGTATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT");
    BOOST_CHECK_EQUAL(read.Qualities, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
    BOOST_CHECK(!read.IsMapped());
    BOOST_CHECK(read.IsFailedQC());
}
