# Written by Scott Norton
# 22-Dec-2022

# Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
# This file is part of DRAGoN.
#
# This source code is licensed under the MIT License found in the
# LICENSE file in the root directory of this source tree.

cmake_minimum_required(VERSION 3.15)
project(drugseq CXX)
set(CMAKE_CXX_STANDARD 23)
include(CTest)

set(DEMULTIPLEX_CXX_SOURCES DemultiplexWorker.cpp FastqReader.cpp BarcodeMatcher.cpp QCStats.cpp QualityFilter.cpp)
add_executable(demultiplex main.cpp ${DEMULTIPLEX_CXX_SOURCES})
target_link_libraries(demultiplex PUBLIC nlohmann_json::nlohmann_json drugseq::drugseq kseq++::kseq++ ${BAMTOOLS_LDFLAGS} Boost::log)
target_compile_options(demultiplex PUBLIC ${BAMTOOLS_CFLAGS})
target_link_options(demultiplex PUBLIC)
install(TARGETS demultiplex)

add_executable(test-demultiplex test_main.cpp ${DEMULTIPLEX_CXX_SOURCES})
target_link_libraries(test-demultiplex PUBLIC drugseq::drugseq kseq++::kseq++ ${BAMTOOLS_LDFLAGS} Boost::unit_test_framework Boost::log)
target_compile_definitions(test-demultiplex PUBLIC BOOST_TEST_MAIN BOOST_TEST_DYN_LINK)
target_compile_options(test-demultiplex PUBLIC ${BAMTOOLS_CFLAGS})
target_link_options(test-demultiplex PUBLIC)
add_test(NAME demultiplex COMMAND test-demultiplex)
