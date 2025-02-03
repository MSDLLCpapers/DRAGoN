#!/usr/bin/env python

"""
Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
This file is part of DRAGoN.

This source code is licensed under the MIT License found in the
LICENSE file in the root directory of this source tree.
"""

import argparse
import json
import os
import re
from collections.abc import Iterable
from typing import TextIO, TypedDict

import pandas as pd
import scipy.io as sio
import scipy.sparse as sp


class DemuxReport(TypedDict):
    read_count: int
    well_counts: dict[str, int]
    written_counts: list[int]
    mismatch_counts: list[dict[str, int]]
    fail_counts: dict[str, dict[str, int]]
    unmatched_counts: dict[str, int]


def load_trimming_log(filenames: list[str]):
    index = pd.Index(
        [
            "Total reads processed",
            "Reads with adapters",
            "Reads written (passing filters)",
            "Total basepairs processed",
            "Quality-trimmed",
            "Total written (filtered)",
        ]
    )
    if filenames is None:
        return pd.Series(0, index=index)

    def inner(filename: str) -> pd.Series:
        ret = pd.Series(0, index=index)
        pattern = re.compile("^(?P<key>.+):\s+(?P<count>[\d,]+)")
        with open(filename) as fp:
            for line in fp:
                if m := pattern.match(line):
                    if m["key"] in ret.index:
                        ret[m["key"]] += int(m["count"].replace(",", ""))
        return ret

    df = sum(inner(filename) for filename in filenames)
    df.to_csv("trimming_summary.txt", sep="\t")
    return df


def read_dedup_report(filenames: Iterable[str]) -> pd.DataFrame:
    dedup_stats = sum(pd.read_csv(fname, sep="\t", index_col=0) for fname in filenames)
    dedup_stats.index.name = "barcode"
    dedup_stats.to_csv("dedup_stats.tsv", sep="\t")
    return dedup_stats


def read_featurecounts(filenames: list[str]) -> pd.DataFrame:
    fcntsum = sum(
        pd.read_csv(fname, sep="\t", index_col=0, header=0, names=["Status", "Count"])
        for fname in filenames
    )
    fcntsum.to_csv("featureCounts.summary", sep="\t")
    return fcntsum


def read_counts_matrices(
    dirnames: list[str], matrix_names: list[str]
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    bcs = pd.read_csv(f"{dirnames[0]}/barcodes.tsv", sep="\t", header=None)
    fts = pd.read_csv(
        f"{dirnames[0]}/features.tsv",
        sep="\t",
        header=None,
        names=["feature", "name", "kind"],
    )
    os.makedirs("DRAGoN.out", exist_ok=True)
    bcs.to_csv("DRAGoN.out/barcodes.tsv", header=False, index=False, sep="\t")
    fts.to_csv("DRAGoN.out/features.tsv", header=False, index=False, sep="\t")
    dedup_report = read_dedup_report((f"{dname}/dedup_stats.tsv" for dname in dirnames))

    def read_sparse_matrix(fname):
        try:
            return sio.mmread(fname).todense()
        except Exception as e:
            raise IOError(fname + " is malformatted") from e

    def inner(name):
        matrix_name = "matrix" if name == "Unique" else f"UniqueAndMult-{name}"
        mtx = sum(
            read_sparse_matrix(f"{dirname}/{matrix_name}.mtx") for dirname in dirnames
        )
        sio.mmwrite(f"DRAGoN.out/{matrix_name}.mtx", sp.csc_matrix(mtx))
        frame = pd.DataFrame(
            mtx, fts.feature, dedup_report.index.drop("noBC", errors="ignore")
        )
        counts_name = matrix_name.replace("UniqueAndMult", "").replace("matrix", "")
        frame.to_csv(f"counts{counts_name}.tsv", sep="\t")
        return frame

    return {name: inner(name) for name in matrix_names}, dedup_report


def read_downsampling_report(downsampling_reports_in: list[str]):
    return (
        pd.concat(
            [
                pd.read_csv(
                    fname,
                    sep="\t",
                )
                for fname in downsampling_reports_in
            ],
            axis=0,
        )
        .groupby(["sample", "index", "well_name"], as_index=False)
        .agg(
            replicate=pd.NamedAgg("replicate", "|".join),
            qc_pass_reads=pd.NamedAgg("qc_pass_reads", "sum"),
            limit=pd.NamedAgg("limit", lambda x: "|".join(map(str, x))),
            seed=pd.NamedAgg("seed", lambda x: "|".join(map(str, x))),
            fraction=pd.NamedAgg("fraction", lambda x: "|".join(map(str, x))),
            out=pd.NamedAgg("out", "sum"),
        )
        .drop(["sample", "index", "replicate"], axis=1)
        .set_index("well_name")
    )


def read_unmapped_report(unmapped_report_in: str):
    return pd.read_csv(
        unmapped_report_in,
        sep="\t",
        index_col=0,
        names=["bcid", "count"],
        dtype={"bcid": str},
    )


def compile_report(
    demux: DemuxReport,
    dedup: pd.DataFrame,
    unmapped_report: pd.DataFrame,
    downsample_report: pd.DataFrame,
):
    """
    Consolidate components of the output statistics into a single table.
    :param demux: JSON-serializable dict from demux.json
    :param dedup: DataFrame with deduplication QC stats
    :return: DataFrame combining the demux, dedup, and counting columns
    """
    mismatch = len(demux["mismatch_counts"]) - 1
    total_in = pd.Series(demux["well_counts"], name="TOTAL_IN")
    fail_mat = pd.DataFrame(demux["fail_counts"])
    mismatch_mat = pd.DataFrame(
        demux["mismatch_counts"],
        index=["Exact"] + [f"{i + 1}MM" for i in range(mismatch)],
    ).transpose()
    unmapped_report = unmapped_report.reindex(
        [f"{i:04d}" for i in range(len(mismatch_mat))] + ["noBC"], fill_value=0
    )
    unmapped_report["index"] = total_in.index
    unmapped_report = unmapped_report.set_index("index")
    if (dedup["DISCARD_UNMAPPED"] == 0).all():
        dedup["DISCARD_UNMAPPED"] = unmapped_report["count"].rename("DISCARD_UNMAPPED")
    filter_remaining = pd.DataFrame(
        {
            "PCT_LOSS_QC": total_in,
            "PCT_LOSS_DOWNSAMPLING": fail_mat["PASS"],
            "PCT_LOSS_UNMAPPED": downsample_report["out"],
            "PCT_LOSS_UNASSIGNED": downsample_report["out"] - unmapped_report["count"],
            "PCT_LOSS_UMI": dedup["PASSED_QC"],
            "FINAL_UMI": dedup["UMI_COUNTS"],
        }
    )
    ret = pd.concat(
        [
            total_in,
            (filter_remaining.diff(periods=-1, axis=1) * 100 / filter_remaining).drop(
                "FINAL_UMI", axis=1
            ),
            dedup["UMI_COUNTS"].rename("FINAL_UMI"),
            fail_mat.rename(columns=lambda x: f"QC_{x.upper()}"),
            mismatch_mat.rename(columns=lambda x: f"MATCH_{x.upper()}"),
            downsample_report.rename(columns=lambda x: f"DOWNSAMPLE_{x.upper()}"),
            unmapped_report.rename(columns=lambda x: f"UNMAP_{x.upper()}"),
            dedup.rename(columns=lambda x: f"DEDUP_{x.upper()}"),
        ],
        axis=1,
    )
    ret.index.name = "barcode"
    ret.to_csv("DRAGoNreport.txt", sep="\t")
    return ret


class CLI(argparse.Namespace):
    matrices: list[str]
    dmux: TextIO
    features: list[str]
    ambiguous: list[str]
    downsampling_report: list[str]
    unmapped_report: str

    def __init__(self, args=None):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--matrices",
            required=True,
            nargs="+",
            help="Directories containing output matrices",
        )
        parser.add_argument(
            "--dmux",
            required=True,
            type=argparse.FileType(),
            help="JSON report from merge_demux_json.py (Preprocessing:MergeDemuxJSON)",
        )
        parser.add_argument(
            "--features",
            required=True,
            nargs="+",
            help="TSV reports from featureCounts (DRAGoN:FeatureCounts)",
        )
        parser.add_argument(
            "--ambiguous",
            required=True,
            nargs="+",
            help="Strategies used to distribute UMIs assigned to two or more genes",
        )
        parser.add_argument(
            "--downsample-report",
            dest="downsampling_report",
            required=True,
            nargs="+",
            help="Result from downsampling reads",
        )
        parser.add_argument(
            "--unmapped-report",
            dest="unmapped_report",
            required=True,
            help="Result from collating unmapped reads",
        )
        _, self.remaining_args = parser.parse_known_args(args, self)

    def main(self):
        demux: DemuxReport = json.load(self.dmux)
        read_featurecounts(self.features)
        _, dedup = read_counts_matrices(self.matrices, self.ambiguous)
        downsampling_report = read_downsampling_report(self.downsampling_report)
        unmapped_report = read_unmapped_report(self.unmapped_report)
        compile_report(demux, dedup, unmapped_report, downsampling_report)


if __name__ == "__main__":
    CLI().main()
