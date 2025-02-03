"""
Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
This file is part of DRAGoN.

This source code is licensed under the MIT License found in the
LICENSE file in the root directory of this source tree.
"""

import os

__version__ = "2025-01-03"
try:
    with open(f"{os.path.dirname(__file__)}/.commit") as fp:
        __version__ = fp.read().strip() or __version__
except OSError:
    __version__ = os.getenv("NF_WORKFLOW_VERSION", __version__)
