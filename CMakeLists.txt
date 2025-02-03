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
message(STATUS CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH})

include(CTest)

find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(PkgConfig REQUIRED)
cmake_policy(PUSH)
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.30)
cmake_policy(VERSION 3.30)
cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost REQUIRED COMPONENTS unit_test_framework log)
cmake_policy(POP)
find_package(kseq++ REQUIRED)
find_package(Eigen3 REQUIRED MODULE)
pkg_check_modules(BAMTOOLS REQUIRED bamtools-1)

add_subdirectory(src/libdrugseq)

set(PROJECT_EXECUTABLES "demultiplex;deduplicate")
foreach(EXECUTABLE IN LISTS PROJECT_EXECUTABLES)
    add_subdirectory("src/${EXECUTABLE}")
    add_dependencies(${EXECUTABLE} drugseq)
endforeach()
