// Written by Scott Norton
// 14-Feb-2023

// Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
// This file is part of DRAGoN.
//
// This source code is licensed under the MIT License found in the
// LICENSE file in the root directory of this source tree.

#include <boost/test/unit_test.hpp>
#include <drugseq.hpp>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace drugseq;

// Unit tests start here
BOOST_AUTO_TEST_CASE(test_strjoin) {
    BOOST_CHECK_EQUAL("hello world", strjoin({"hello", "world"}, " "));
}

BOOST_AUTO_TEST_CASE(test_strsplit) {
    std::vector<std::string> result = strsplit("hello: world :", ": ");
    std::vector<std::string> expected {"hello", "world :"};
    BOOST_CHECK_EQUAL_COLLECTIONS(result.begin(), result.end(), expected.begin(), expected.end());
}

BOOST_AUTO_TEST_CASE(test_shlexjoin) {
    BOOST_CHECK_EQUAL("hello \"world\"", shlexjoin({"hello", "\"world\""}));
    BOOST_CHECK_EQUAL("hello \"world foo\"", shlexjoin({"hello", "world foo"}));
    BOOST_CHECK_EQUAL("hello \"world \\\"foo\\\"\"", shlexjoin({"hello", "world \"foo\""}));
}
