#!/usr/bin/env bash

HOME_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

set -x
rm -rf "${HOME_DIR}/bin/demultiplex" "${HOME_DIR}/bin/deduplicate" "${HOME_DIR}/build"
