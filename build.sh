#!/usr/bin/env bash

HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

set -xe
mkdir -p build
cmake "$HOME_DIR" -B build -DCMAKE_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX:-$HOME_DIR}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="${CONDA_ROOT}" -DCMAKE_MODULE_PATH=${CONDA_PREFIX}/unpacked_source/cmake -DDRUGSEQ_VERSION_STR="$(cat /version.txt)"
cmake --build build -j 4
cmake --install build --strip
