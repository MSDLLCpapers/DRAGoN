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

set(DEDUPLICATE_CXX_SOURCES DeduplicateWorker.cpp DuplicateDetector.cpp SamReader.cpp)
add_executable(deduplicate main.cpp ${DEDUPLICATE_CXX_SOURCES})
target_link_libraries(deduplicate PUBLIC Threads::Threads drugseq::drugseq ${BAMTOOLS_LDFLAGS})
target_include_directories(deduplicate PUBLIC ${EIGEN3_INCLUDE_DIR})
target_compile_options(deduplicate PUBLIC ${BAMTOOLS_CFLAGS})
target_link_options(deduplicate PUBLIC)
install(TARGETS deduplicate)

add_executable(test-deduplicate test_main.cpp ${DEDUPLICATE_CXX_SOURCES})
target_link_libraries(test-deduplicate PUBLIC drugseq::drugseq Boost::unit_test_framework Boost::log Threads::Threads ${BAMTOOLS_LDFLAGS})
target_include_directories(test-deduplicate PUBLIC ${EIGEN3_INCLUDE_DIR})
target_compile_definitions(test-deduplicate PUBLIC BOOST_TEST_MAIN BOOST_TEST_DYN_LINK)
target_compile_options(test-deduplicate PUBLIC ${BAMTOOLS_CFLAGS})
target_link_options(test-deduplicate PUBLIC)
add_test(NAME deduplicate COMMAND test-deduplicate)
