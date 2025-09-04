#!/usr/bin/env python

"""
Written by Scott Norton
2025-02-27

Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
This file is part of TDP-DISCOVERY-INFORMATICS-DRUG-SEQ-PIPELINE.

This source code is licensed under the MIT License found in the
LICENSE file in the root directory of this source tree.
"""

import argparse
import csv
import json
import typing


def quantile(nums, quantile, is_sorted=True):
    if not is_sorted:
        nums = sorted(nums)
    index_f = quantile * len(nums)
    index = int(index_f)
    if index == index_f:
        return nums[index]
    else:
        frac = index_f - index
        return nums[index] * (1 - frac) + nums[index + 1] * frac


class CLI(argparse.Namespace):
    barcodes: typing.TextIO
    demux_log: typing.TextIO

    def __init__(self, args=None):
        parser = argparse.ArgumentParser()
        parser.add_argument("barcodes", type=argparse.FileType())
        parser.add_argument("demux_log", type=argparse.FileType())
        parser.parse_args(args, self)

    def main(self):
        json_data = json.load(self.demux_log)
        total_size = json_data["read_count"]
        if total_size == 0:
            raise ValueError("no reads in experiment")
        well_sizes = sorted(
            json_data["fail_counts"]["PASS"][row[1]]
            for row in csv.reader(self.barcodes, dialect="excel-tab")
            if len(row) < 5 or row[4] != "True"
        )
        if sum(well_sizes) == 0:
            raise ValueError("no QC pass reads in unmasked wells")
        q25 = quantile(well_sizes, 0.25)
        q75 = quantile(well_sizes, 0.75)
        thresh = q75 + 1.75 * (q75 - q25)
        return min(thresh, 0.15 * total_size)


if __name__ == "__main__":
    cli = CLI()
    result = cli.main()
    print(result)
