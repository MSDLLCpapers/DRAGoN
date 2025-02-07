#!/usr/bin/env python

"""
Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
This file is part of DRAGoN.

This source code is licensed under the MIT License found in the
LICENSE file in the root directory of this source tree.
"""

import abc
import argparse
import collections
import functools
import pathlib
import re
from typing import Literal, TextIO

import pandas as pd
from numpy import inf
from scipy.special import logsumexp

M_LN1BIL = 20.72326583694641  # log(10e9)
M_LN1MIL = 13.815510557964274  # log(10e6)


@functools.cache
def gene_effective_length_from_gtf(
    pth: pathlib.Path, attr_key="gene_id", exon_feature="exon"
):
    exons = collections.defaultdict(set)
    with open(pth) as fp:
        for line in fp:
            if line.startswith("#"):
                continue
            seqname, source, feature, start, end, score, strand, frame, attribute = (
                line.split("\t")
            )
            if feature == exon_feature:
                gene_id = re.search(f'{attr_key} "(.+?)";', attribute)[1]
                exons[gene_id].add((int(start), int(end)))

    # adapted from https://stackoverflow.com/a/59401068
    # Collapse overlapping exons, this could arise i.e. from an alternative isoform with an A3'SS, A5'SS, or IR event
    return pd.Series(
        {
            name: sum(
                e - s
                for s, e in functools.reduce(
                    lambda acc, el: (
                        acc[:-1:] + [(min(*acc[-1], *el), max(*acc[-1], *el))]
                        if acc[-1][1] >= el[0] - 1
                        else acc + [el]
                    ),
                    sc[1:],
                    sc[:1],
                )
            )
            for name, sc in ((name, sorted(coords)) for name, coords in exons.items())
        }
    )


class CountsNormalizer(abc.ABC):
    def __init__(
        self, gtffile: pathlib.Path, attribute: str = "gene_id", feature: str = "exon"
    ):
        self.effective_lengths = gene_effective_length_from_gtf(
            gtffile, attr_key=attribute, exon_feature=feature
        )

    @abc.abstractmethod
    def adjustment_impl(
        self, raw_counts: pd.DataFrame, effective_lengths: pd.DataFrame
    ) -> pd.DataFrame:
        return NotImplemented

    def __call__(
        self,
        raw_counts: pd.DataFrame,
        fragment_length: list[float] = None,
        fragment_length_names: list[str] = None,
    ) -> pd.DataFrame:
        if fragment_length and fragment_length_names:
            effective_lengths = pd.DataFrame(
                {
                    name: self.effective_lengths - fgment_len
                    for name, fgment_len in zip(fragment_length_names, fragment_length)
                }
            ).clip(0)
        else:
            effective_lengths = pd.DataFrame(
                {name: self.effective_lengths for name in raw_counts}
            )
        result = self.adjustment_impl(raw_counts, effective_lengths)
        result[effective_lengths <= 0] = 0
        return result.fillna(0).reindex(index=self.effective_lengths.index)


class TPMNormalizer(CountsNormalizer):
    def adjustment_impl(
        self, raw_counts: pd.DataFrame, effective_lengths: pd.DataFrame
    ) -> pd.DataFrame:
        rate = (
            raw_counts.transform("log") - effective_lengths.transform("log")
        ).fillna(-inf)
        denom = rate.apply(logsumexp)
        return (rate - denom + M_LN1MIL).transform("exp")


class FPKMNormalizer(CountsNormalizer):
    def adjustment_impl(
        self, raw_counts: pd.DataFrame, effective_lengths: pd.DataFrame
    ) -> pd.DataFrame:
        rate = raw_counts.transform("log").fillna(-inf)
        denom = rate.apply(logsumexp)
        return (rate - denom - effective_lengths.transform("log") + M_LN1BIL).transform(
            "exp"
        )


class CLI(argparse.Namespace):
    mode: Literal["tpm", "fpkm"]
    raw_counts: pd.DataFrame
    _gtf: pathlib.Path
    fragment_length: list[float]
    fragment_length_names: list[str]
    outfile: TextIO
    attribute: str
    exon_feature: str

    _parser = argparse.ArgumentParser()
    _parser.add_argument("mode", choices=["tpm", "fpkm"])
    _parser.add_argument(
        "-c",
        dest="raw_counts",
        type=lambda x: pd.read_csv(x, sep="\t", index_col=0),
        required=True,
    )
    _parser.add_argument("-a", dest="_gtf", type=pathlib.Path, required=True)
    _parser.add_argument("-f", dest="fragment_length", type=float, nargs="+")
    _parser.add_argument("-n", dest="fragment_length_names", nargs="+")
    _parser.add_argument(
        "-o", dest="outfile", type=argparse.FileType("w"), required=True
    )
    _parser.add_argument("-i", dest="attribute", default="gene_id")
    _parser.add_argument("-e", dest="exon_feature", default="exon")

    def __init__(self, args=None):
        self.__class__._parser.parse_args(args, self)

    def main(self):
        cls: type[CountsNormalizer] = {"fpkm": FPKMNormalizer, "tpm": TPMNormalizer}[
            self.mode
        ]
        cls(self._gtf, self.attribute, self.exon_feature)(
            self.raw_counts, self.fragment_length, self.fragment_length_names
        ).to_csv(
            self.outfile,
            sep="\t",
            float_format="%.10g",  # workaround for GitHub Actions precision issues
        )


if __name__ == "__main__":
    CLI().main()
