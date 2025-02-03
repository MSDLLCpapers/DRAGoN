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

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/drugseq/__version.hpp")
    MESSAGE(FATAL_ERROR "File \"${CMAKE_CURRENT_SOURCE_DIR}/include/drugseq/__version.hpp\" found. Its presence breaks the build system. Please remove at once.")
endif()

execute_process(
    COMMAND git show -s --format="%cs+g%h"
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    OUTPUT_VARIABLE DRUGSEQ_VERSION_STR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_custom_target(version_header DEPENDS "include/drugseq/__version.hpp")
add_custom_command(COMMAND touch "${CMAKE_CURRENT_SOURCE_DIR}/include/drugseq/version.hpp" OUTPUT "include/drugseq/version.hpp" "include/drugseq/__version.hpp")

set(LIBDRUGSEQ_CXX_SOURCES Gtf.cpp include/drugseq/version.hpp)
add_library(drugseq ${LIBDRUGSEQ_CXX_SOURCES})
add_library(drugseq::drugseq ALIAS drugseq)
target_link_libraries(drugseq PUBLIC ZLIB::ZLIB ${BAMTOOLS_LDFLAGS})
target_compile_definitions(drugseq PUBLIC DRUGSEQ_VERSION_STR="\"${DRUGSEQ_VERSION_STR}\"")
target_compile_options(drugseq PUBLIC ${BAMTOOLS_CFLAGS})
target_include_directories(drugseq PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_features(drugseq PUBLIC cxx_std_17)
target_link_options(drugseq PUBLIC)
set(DRUGSEQ_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include;${ZSTR_INCLUDE_DIRS}" CACHE PATH "" FORCE)
set(DRUGSEQ_LIBRARIES "${ZSTR_LIBRARIES};${BAMTOOLS_LIBRARIES}" CACHE PATH "" FORCE)

add_executable(test-libdrugseq test_main.cpp ${LIBDRUGSEQ_CXX_SOURCES})
target_link_libraries(test-libdrugseq PUBLIC drugseq::drugseq Boost::unit_test_framework)
target_compile_definitions(test-libdrugseq PUBLIC BOOST_TEST_MAIN BOOST_TEST_DYN_LINK)
add_test(NAME libdrugseq COMMAND test-libdrugseq)
