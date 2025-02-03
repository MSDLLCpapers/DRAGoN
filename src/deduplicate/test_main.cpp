// Written by Scott Norton
// 15-Feb-2023

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#include <boost/test/unit_test.hpp>
#include <boost/log/trivial.hpp>
#include <filesystem>
#include <vector>
#include <string>
#include <unordered_set>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/SparseExtra>
#include "DeduplicateWorker.hpp"

namespace fs = std::filesystem;

fs::path __file__ {__FILE__};
fs::path data_dir = __file__.parent_path().parent_path() / "test-data";
fs::path inbam = data_dir / "tmp1.bam",
            refbam = data_dir / "tmp1_MarkDups_ref.bam",
            gtf = data_dir / "genes.gtf",
            barcodes = data_dir / "Simulated_1_barcodes.txt",
            refouts = data_dir / "DRAGoN.out.ref",
            cmpouts = data_dir / "DRAGoN.out";

BOOST_AUTO_TEST_CASE(deduplicate_2000) {
    BOOST_LOG_TRIVIAL(debug) << "data_dir = " << data_dir;
    BOOST_TEST_REQUIRE(fs::exists(data_dir));
    BOOST_TEST_REQUIRE(fs::exists(inbam));
    BOOST_TEST_REQUIRE(fs::exists(refbam));
    BOOST_TEST_REQUIRE(fs::exists(gtf));
    BOOST_TEST_REQUIRE(fs::exists(barcodes));
    BOOST_TEST_REQUIRE(fs::exists(refouts));

    fs::remove_all(cmpouts);
    fs::remove(data_dir / "tmp1_MarkDups.bam");

    std::vector<std::string> cli {
        "deduplicate", inbam, barcodes, gtf, cmpouts,
        "--mmap-strategies", "Unique", "Uniform", "PropUnique", "EM", "Rescue",
        "--keep-bam",
        "--umilen", "12"
    };

    {
        DeduplicateWorker worker {
            cli,
            inbam,
            1,
            barcodes,
            gtf,
            {MultimapStrategy::UNIFORM, MultimapStrategy::PROPUNIQUE, MultimapStrategy::EM, MultimapStrategy::RESCUE},
            12,
            2,
            10,
            100,
            .01,
            .01,
            true
        };
        BOOST_REQUIRE_EQUAL(worker.run(), 0);
        worker.save_output(cmpouts);
    }
}

BOOST_AUTO_TEST_CASE(uniques_matrix, * boost::unit_test::depends_on("deduplicate_2000"))
{
    Eigen::SparseMatrix<unsigned long long> refmat, cmpmat;
    Eigen::loadMarket(refmat, (refouts / "matrix.mtx").string());
    Eigen::loadMarket(cmpmat, (cmpouts / "matrix.mtx").string());
    BOOST_CHECK(refmat.isApprox(cmpmat));
}

void test_multimap_matrix_sub(const std::string& s) {
    Eigen::SparseMatrix<double> refmat, cmpmat;
    std::string mtxname = "UniqueAndMult-" + s + ".mtx";
    Eigen::loadMarket(refmat, (refouts / mtxname).string());
    Eigen::loadMarket(cmpmat, (cmpouts / mtxname).string());
    BOOST_CHECK(refmat.isApprox(cmpmat));
}

BOOST_AUTO_TEST_CASE(uniform_matrix, * boost::unit_test::depends_on("deduplicate_2000")) {
    test_multimap_matrix_sub("Uniform");
}

BOOST_AUTO_TEST_CASE(propunique_matrix, * boost::unit_test::depends_on("deduplicate_2000")) {
    test_multimap_matrix_sub("PropUnique");
}

BOOST_AUTO_TEST_CASE(em_matrix, * boost::unit_test::depends_on("deduplicate_2000")) {
    test_multimap_matrix_sub("EM");
}

BOOST_AUTO_TEST_CASE(rescue_matrix, * boost::unit_test::depends_on("deduplicate_2000")) {
    test_multimap_matrix_sub("Rescue");
}
