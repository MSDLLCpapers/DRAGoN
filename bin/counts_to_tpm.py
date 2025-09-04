#!/usr/bin/env python

"""
Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
This file is part of DRAGoN.

This source code is licensed under the MIT License found in the
LICENSE file in the root directory of this source tree.
"""

# Performs TPM and FPKM normalization of expression counts matrices
# Scott Norton
# 2023-05-02
# Version: 2025-01-21

import abc
import argparse
import collections
import functools
import pathlib
from typing import Literal, TextIO

import pandas as pd

try:
    from numpy import inf
except ImportError:
    from numpy import infty as inf

from libstarburst import read_gtf
from scipy.special import logsumexp

M_LN1BIL = 20.72326583694641  # log(10e9)
M_LN1MIL = 13.815510557964274  # log(10e6)


@functools.cache
def gene_effective_length_from_gtf(
    pth: pathlib.Path,
    attr_key="gene_id",
    exon_feature="exon",
    biotype_attr="gene_biotype",
    filter_rrna=False,
):
    exons = collections.defaultdict(set)
    denominator_genes: set[str] = set()
    for _, _, feature, start, end, _, _, _, attribute in read_gtf(pth):
        if attr_key not in attribute:
            continue
        gene_id = attribute[attr_key]
        cur_exons = exons[gene_id]
        if feature == exon_feature:
            cur_exons.add((start, end))
        if not (filter_rrna and "rRNA" in attribute.get(biotype_attr, "")):
            denominator_genes.add(gene_id)

    # adapted from https://stackoverflow.com/a/59401068
    # Collapse overlapping exons, this could arise i.e. from an alternative isoform with an A3'SS, A5'SS, or IR event
    return (
        pd.Series(
            {
                name: (
                    sum(
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
                )
                for name, sc in (
                    (name, sorted(coords)) for name, coords in exons.items()
                )
            }
        ),
        pd.Index(denominator_genes),
    )


class CountsNormalizer(abc.ABC):
    def __init__(
        self,
        gtffile: pathlib.Path,
        attribute: str = "gene_id",
        feature: str = "exon",
        biotype_attr: str = "gene_biotype",
        filter_rrna=False,
    ):
        self.effective_lengths, self.denominator_genes = gene_effective_length_from_gtf(
            gtffile,
            attr_key=attribute,
            exon_feature=feature,
            biotype_attr=biotype_attr,
            filter_rrna=filter_rrna,
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
        result = self.adjustment_impl(
            raw_counts.reindex(index=self.effective_lengths.index, fill_value=0),
            effective_lengths,
        )
        result[effective_lengths <= 0] = 0
        return result.fillna(0)


class TPMNormalizer(CountsNormalizer):
    def adjustment_impl(
        self, raw_counts: pd.DataFrame, effective_lengths: pd.DataFrame
    ) -> pd.DataFrame:
        rate = (
            raw_counts.transform("log") - effective_lengths.transform("log")
        ).fillna(-inf)
        denom = rate.loc[self.denominator_genes].apply(logsumexp)
        return (rate - denom + M_LN1MIL).transform("exp")


class FPKMNormalizer(CountsNormalizer):
    def adjustment_impl(
        self, raw_counts: pd.DataFrame, effective_lengths: pd.DataFrame
    ) -> pd.DataFrame:
        rate = raw_counts.transform("log").fillna(-inf)
        denom = rate.loc[self.denominator_genes].apply(logsumexp)
        return (rate - denom - effective_lengths.transform("log") + M_LN1BIL).transform(
            "exp"
        )


class CLI(argparse.Namespace):
    mode: Literal["tpm", "fpkm"]
    raw_counts: pd.DataFrame
    gtf: pathlib.Path
    fragment_length: list[float]
    fragment_length_names: list[str]
    outfile: TextIO
    attribute: str = "gene_id"
    exon_feature: str = "exon"
    biotype_attr: str = "gene_biotype"
    filter_rrna: bool = False

    def __init__(self, args=None):

        parser = argparse.ArgumentParser()
        parser.add_argument("mode", choices=["tpm", "fpkm"])
        parser.add_argument(
            "-c",
            dest="raw_counts",
            type=lambda x: pd.read_csv(x, sep="\t", index_col=0),
            required=True,
        )
        parser.add_argument("-a", dest="gtf", type=pathlib.Path, required=True)
        parser.add_argument("-f", dest="fragment_length", type=float, nargs="+")
        parser.add_argument("-n", dest="fragment_length_names", nargs="+")
        parser.add_argument(
            "-o", dest="outfile", type=argparse.FileType("w"), required=True
        )
        parser.add_argument("-i", dest="attribute", default=CLI.attribute)
        parser.add_argument("-e", dest="exon_feature", default=CLI.exon_feature)
        parser.add_argument("-b", dest="biotype_attr", default=CLI.biotype_attr)
        parser.add_argument("-R", dest="filter_rrna", action="store_true")
        parser.parse_args(args, self)

    def main(self):
        cls: type[CountsNormalizer] = {"fpkm": FPKMNormalizer, "tpm": TPMNormalizer}[
            self.mode
        ]
        norm = cls(
            self.gtf,
            attribute=self.attribute,
            feature=self.exon_feature,
            biotype_attr=self.biotype_attr,
            filter_rrna=self.filter_rrna,
        )
        norm(
            self.raw_counts,
            fragment_length=self.fragment_length,
            fragment_length_names=self.fragment_length_names,
        ).rename_axis(self.attribute, axis=0).to_csv(
            self.outfile, sep="\t", float_format="%.8g"
        )


if __name__ == "__main__":
    CLI().main()
