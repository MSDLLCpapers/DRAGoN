#!/usr/bin/env python

"""
Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
This file is part of DRAGoN.

This source code is licensed under the MIT License found in the
LICENSE file in the root directory of this source tree.
"""

import argparse
import json
import pathlib
import re
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


def load_dedup_report(flagstat_fn: list[str]) -> tuple[pd.Series, pd.DataFrame]:
    """
    Loads the counts matrices and flagstat reports from STARsolo
    :param dedup_fn: Path to counts matrix
    :param flagstat_fn: List of paths to samtools flagstat jsons
    :return: 2-tuple with a pd.Series of summary counts and a DataFrame combining the QC-passed read counts with flags. Additionally, dumps the summed flagstat jsons to a new file 'flagstat.json'.
    """
    json_orig = [pd.read_json(fn) for fn in flagstat_fn]
    consolidated: pd.DataFrame = sum(json_orig)
    json_qcpass = pd.concat([x["QC-passed reads"] for x in json_orig])
    consolidated.loc["mapped %"] = (
        consolidated.loc["mapped"] * 100 / consolidated.loc["total"]
    )
    consolidated.loc["primary mapped %"] = (
        consolidated.loc["primary mapped"] * 100 / consolidated.loc["primary"]
    )
    consolidated.loc["properly paired %"] = (
        consolidated.loc["properly paired"]
        * 100
        / consolidated.loc["paired in sequencing"]
    )
    consolidated.loc["singletons %"] = (
        consolidated.loc["singletons"] * 100 / consolidated.loc["paired in sequencing"]
    )
    consolidated.to_json("flagstat.json", indent=1)
    return json_qcpass


def read_counts_matrices(
    solo_outs: list[str], matrix_names: list[str], barcodes_tsv: pd.DataFrame
):
    bcs = pd.read_csv(f"{solo_outs[0]}/Gene/raw/barcodes.tsv", sep="\t", header=None)
    fts = pd.read_csv(
        f"{solo_outs[0]}/Gene/raw/features.tsv",
        sep="\t",
        header=None,
        names=["feature", "name", "kind"],
    )
    pathlib.Path("Solo.out/Gene/raw").mkdir(parents=True, exist_ok=True)
    bcs.to_csv("Solo.out/Gene/raw/barcodes.tsv", header=False, index=False, sep="\t")
    fts.to_csv("Solo.out/Gene/raw/features.tsv", header=False, index=False, sep="\t")

    # merge Barcodes.stats, Features.stats, Summary.csv
    with open("Solo.out/Barcodes.stats", "w") as bcs:
        for name, value in sum(
            pd.read_fwf(
                f"{dirname}/Barcodes.stats", widths=[50, 15], index_col=0, header=None
            )
            for dirname in solo_outs
        ).itertuples():
            print(f"{name:>50s}{value:>15d}", file=bcs)
    with open("Solo.out/Gene/Features.stats", "w") as fcs:
        for name, value in sum(
            pd.read_fwf(
                f"{dirname}/Gene/Features.stats",
                widths=[50, 15],
                index_col=0,
                header=None,
            )
            for dirname in solo_outs
        ).itertuples():
            print(f"{name:>50s}{value:>15d}", file=fcs)
    scv = (
        pd.concat(
            [
                pd.read_csv(f"{dirname}/Gene/Summary.csv", index_col=0, header=None)
                for dirname in solo_outs
            ],
            axis=1,
        )
        .replace("NoMulti", 0)
        .astype(float)
    )
    scv = scv.apply(
        lambda row: (
            row.sum()
            if row.name == "Number of Reads"
            else (row * scv.loc["Number of Reads"]).sum()
            / scv.loc["Number of Reads"].sum()
        ),
        axis=1,
    )
    scv.to_csv("Solo.out/Gene/Summary.csv", header=False)

    def read_sparse_matrix(fname):
        try:
            return sio.mmread(fname).todense()
        except Exception as e:
            raise IOError(fname + " is malformatted") from e

    def inner(name):
        matrix_name = "matrix" if name == "Unique" else f"UniqueAndMult-{name}"
        mtx = sum(
            read_sparse_matrix(f"{dirname}/Gene/raw/{matrix_name}.mtx")
            for dirname in solo_outs
        )
        sio.mmwrite(f"Solo.out/Gene/raw/{matrix_name}.mtx", sp.csc_matrix(mtx))
        frame = pd.DataFrame(mtx, fts.feature, barcodes_tsv.index)
        counts_name = matrix_name.replace("UniqueAndMult", "").replace("matrix", "")
        frame.to_csv(f"counts{counts_name}.tsv", sep="\t")
        return frame

    return {name: inner(name) for name in matrix_names}


def read_star_report(star_report_in: list[str]):
    star_report_in.sort()
    pat = re.compile(r"^\w+\.(\w{4})_Log\.final\.out")
    df = pd.DataFrame(
        0,
        index=[pat.match(name)[1] for name in star_report_in],
        columns=["STAR_READS_IN", "UNIQUELY_MAPPED", "MULTIMAPPED"],
    )
    for filename in star_report_in:
        index = pat.match(filename)[1]
        with open(filename) as fp:
            for line in fp:
                if "Number of input reads" in line:
                    df.loc[index, "STAR_READS_IN"] += int(
                        line.rstrip("\n").rsplit(None, 1)[1]
                    )
                elif "Uniquely mapped reads number" in line:
                    df.loc[index, "UNIQUELY_MAPPED"] += int(
                        line.rstrip("\n").rsplit(None, 1)[1]
                    )
                elif "Number of reads mapped to multiple loci" in line:
                    df.loc[index, "MULTIMAPPED"] += int(
                        line.rstrip("\n").rsplit(None, 1)[1]
                    )
    return df


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
        unmapped_report_in, sep="\t", index_col=0, names=["index", "count"]
    )


def compile_report(
    demux: DemuxReport,
    frames: dict[str, pd.DataFrame],
    unmapped_report: pd.DataFrame,
    downsample_report: pd.DataFrame,
):
    """
    Consolidate components of the output statistics into a single table.
    :param demux: JSON-serializable dict from demux.json
    :param frames: Final counts matrices as a dict of DataFrame
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
    solo_counts = pd.DataFrame(
        {name: frame.sum(axis=0).round() for name, frame in frames.items()}
    )
    filter_remaining = pd.DataFrame(
        {
            "PCT_LOSS_QC": total_in,
            "PCT_LOSS_DOWNSAMPLING": fail_mat["PASS"],
            "PCT_LOSS_UNMAPPED": downsample_report["out"],
            "PCT_LOSS_UNASSIGNED_UMI": downsample_report["out"]
            - unmapped_report["count"],
            "FINAL_UMI": solo_counts.iloc[:, -1],
        }
    )
    ret = pd.concat(
        [
            total_in,
            (filter_remaining.diff(periods=-1, axis=1) * 100 / filter_remaining).drop(
                "FINAL_UMI", axis=1
            ),
            solo_counts.iloc[:, -1].rename("FINAL_UMI"),
            fail_mat.rename(columns=lambda x: f"QC_{x.upper()}"),
            mismatch_mat.rename(columns=lambda x: f"MATCH_{x.upper()}"),
            downsample_report.rename(columns=lambda x: f"DOWNSAMPLE_{x.upper()}"),
            unmapped_report.rename(columns=lambda x: f"UNMAP_{x.upper()}"),
            solo_counts.rename(columns=lambda x: f"SOLO_{x.upper()}"),
        ],
        axis=1,
    )
    ret.index.name = "barcode"
    ret.to_csv("STARsoloReport.txt", sep="\t")
    return ret


class CLI(argparse.Namespace):
    demux: TextIO
    solo_outs: list[str]
    barcodes: str
    ambiguous: list[str]
    downsampling_report: list[str]
    unmapped_report: str

    def __init__(self, args=None):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--demux",
            required=True,
            type=argparse.FileType(),
            help="JSON report from merge_demux_json.py (Preprocessing:MergeDemuxJSON)",
        )
        parser.add_argument(
            "--solo-outs",
            required=True,
            nargs="+",
            help="Directories containing output matrices",
        )
        parser.add_argument(
            "--barcodes", required=True, help="File containing barcodes definitions"
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
        barcodes_tsv = pd.read_csv(
            self.barcodes,
            sep="\t",
            header=None,
            index_col=1,
            names=["barcode", "well", "x", "y"],
        )
        demux_json: DemuxReport = json.load(self.demux)
        frames = read_counts_matrices(self.solo_outs, self.ambiguous, barcodes_tsv)
        # dedup_json = load_dedup_report(args.flagstat)
        downsampling_report = read_downsampling_report(self.downsampling_report)
        unmapped_report = read_unmapped_report(self.unmapped_report)
        compile_report(demux_json, frames, unmapped_report, downsampling_report)


if __name__ == "__main__":
    CLI().main()
