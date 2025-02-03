#!/usr/bin/env python

"""
Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
This file is part of DRAGoN.

This source code is licensed under the MIT License found in the
LICENSE file in the root directory of this source tree.
"""

import argparse
import array
import collections
import contextlib
import dataclasses
import enum
import json
import logging
import os.path
import shlex
import sys
import typing
from collections.abc import Callable, Iterator
from functools import cache, cached_property, lru_cache

import pandas as pd
import pysam
import regex
from cutadapt.adapters import BackAdapter
from version import __version__


# Dataclass for tracking QC statistics
@dataclasses.dataclass(slots=True)
class ReadCounts:
    read_count: int
    well_counts: dict[str, int]
    written_counts: list[int]
    mismatch_counts: list[dict[str, int]]
    fail_counts: dict[str, dict[str, int]]
    unmatched_counts: collections.Counter[str]


class QCFailReason(enum.Enum):
    PASS = "pass"
    FAIL_UMI_QUALITY = "fail_umi_quality"
    FAIL_AT_CONTENT = "fail_at_content"
    FAIL_GC_CONTENT = "fail_gc_content"
    FAIL_LENGTH = "fail_length"
    FAIL_HOMOPOLYMER = "fail_homopolymer"
    FAIL_ADAPTER_TRIMMING = "fail_adapter_trimming"
    FAIL_QUALITY_TRIMMING = "fail_quality_trimming"
    FAIL_DUPLICATE = "fail_duplicate"


@dataclasses.dataclass(slots=True, unsafe_hash=True)
class ReadPair:
    name: str = dataclasses.field(hash=False)
    cdna_seq: str
    cdna_qual: array.array = dataclasses.field(hash=False)
    bcseq: str
    bcqual: str = dataclasses.field(hash=False)
    umiseq: str
    umiqual: str = dataclasses.field(hash=False)
    bcqual_a: array.array = dataclasses.field(hash=False)
    umiqual_a: array.array = dataclasses.field(hash=False)

    def trim(self, end: int):
        self.cdna_seq = self.cdna_seq[:end]
        self.cdna_qual = self.cdna_qual[:end]


class SamWriter:
    __slots__ = ("outbam", "record", "num_written")

    def __init__(
        self,
        filename: str,
        header: typing.Union[dict, pysam.AlignmentHeader, None],
        cpus: int = 1,
    ):
        self.outbam = pysam.AlignmentFile(filename, "wb", header=header, threads=cpus)
        self.record = pysam.AlignedSegment()
        self.record.is_unmapped = True
        self.num_written = 0

    def write(
        self,
        rp: ReadPair,
        qcpass: QCFailReason,
        bcseq: str,
        mismatch: int,
        fileidx: int,
    ):
        """
        Writes a ReadPair to the appropriate SAM file.
        :param rp: The ReadPair to write
        :param qcpass: Used to determine whether the QCFAIL flag should be set
        :param bcseq: Nucleotide sequence of the true barcode
        :param mismatch: Number of substitutions between the whitelist barcode and the sequenced barcode
        :param fileidx: Integer index of the SAM file to write, specify -1 for "noBC".
        """
        self.record.is_qcfail = qcpass is not QCFailReason.PASS
        self.record.is_duplicate = qcpass is QCFailReason.FAIL_DUPLICATE
        self.record.query_name = rp.name
        self.record.query_sequence = rp.cdna_seq
        self.record.query_qualities = rp.cdna_qual
        rg = f"{fileidx:04d}" if fileidx >= 0 else "noBC"
        self.record.set_tags(
            [
                ("CR", rp.bcseq),
                ("UR", rp.umiseq),
                ("CY", rp.bcqual),
                ("UY", rp.umiqual),
                ("CB", bcseq),
                ("UB", rp.umiseq),
                ("BI", rg),
                ("BM", mismatch),
                ("QC", qcpass.name),
                ("RG", rg),
            ]
        )
        self.outbam.write(self.record)
        self.num_written += 1

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.outbam.__exit__(exc_type, exc_val, exc_tb)


class Demultiplexer:
    IUPAC = "=ACMGRSVTWYHKDBN"
    IUPAC_ARR = (
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        1,
        2,
        4,
        8,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        0,
        15,
        15,
        15,
        1,
        14,
        2,
        13,
        15,
        15,
        4,
        11,
        15,
        15,
        12,
        15,
        3,
        15,
        15,
        15,
        15,
        5,
        6,
        8,
        15,
        7,
        9,
        15,
        10,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        1,
        14,
        2,
        13,
        15,
        15,
        4,
        11,
        15,
        15,
        12,
        15,
        3,
        15,
        15,
        15,
        15,
        5,
        6,
        8,
        15,
        7,
        9,
        15,
        10,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
        15,
    )

    def __init__(
        self,
        barcodes_fn,
        bclen=16,
        umistart=17,
        umilen=10,
        split=0,
        mismatch=2,
        check=True,
        cpus=1,
        min_length=20,
        min_qual=20,
        min_qual_bases=6,
        max_AT=0.9,
        max_GC=0.9,
        min_homop_len=10,
        max_homop_mmatch=0,
        max_bc_ns=2,
        BAM_WRITE_MULTIPLE=False,
        adapter: str = None,
        write_qcfail=True,
    ):
        """Initialize a Demultiplexer instance.

        :param barcodes_fn: Path to a TSV file containing whitelisted barcode sequences and well names
        :param bclen: Length of the barcode. While this can technically be inferred from the barcode filename,
            this implementation does not do that.
        :param umistart: Start position of the UMI sequence within the R1 read, usually one past the barcode end.
        :param umilen: Number of positions in the UMI, usually len(R1) - bclen.
        :param split: Divide barcodes into this many files. Barcodes are distributed evenly, but this may not translate
            to an equal distribution of reads. Set this to 0 to output one file per barcode. Also outputs a no-barcode
            file.
        :param mismatch: Maximum number of mismatches between the sequenced barcode and the whitelist.
        :param check: If True (the default), checks the barcode whitelist for possible collisions, and emits a warning
            if any two distinct barcodes are too similar to each other based on the mismatches parameter.
        :param cpus: Number of compression threads for BAM writing
        :param min_length: Minimum number of cDNA bases after trimming. Default: 20
        :param min_qual: Quality threshold for cDNA and UMI. Default: 20
        :param min_qual_bases: Minimum number of UMI bases that must exceed min_qual. Default: 6
        :param max_AT: Float maximum AT fraction in cDNA. Default: 0.9
        :param max_GC: Float maximum GC fraction in cDNA. Ddefault: 0.9
        :param min_homop_len: Minimum homopolymer to trim from 3' end (default: 10)
        :param max_homop_mmatch: Maximum number of mismatches allowed in homopolymer sequence (default: 0)
        :param max_bc_ns: Maximum number of N bases allowed in a barcode. Any barcode exceeding this limit will be
            treated as "noBC". (default: 2)
        :param BAM_WRITE_MULTIPLE: Output multiple BAM files. When selected, the cpus parameter is ignored.
        :param write_qcfail: Controls whether to write reads that failed QC to the BAM file.
        """

        # Set params
        self.bclen = bclen
        self.umistart = umistart
        self.umilen = umilen
        self.split = split
        self.mismatch = mismatch
        self.min_length = min_length
        self.min_qual = min_qual
        self.min_qual_bases = min_qual_bases
        self.max_AT = max_AT
        self.max_GC = max_GC
        self.min_homop_len = min_homop_len
        self.max_homop_mmatch = max_homop_mmatch
        self.max_bc_ns = max_bc_ns
        self.cpus = cpus
        self.BAM_WRITE_MULTIPLE = BAM_WRITE_MULTIPLE
        self._adapter = adapter
        self.write_qcfail = write_qcfail

        self.init_filters()

        # Read the barcode sequences and well names
        self.barcodes = pd.read_csv(
            barcodes_fn, sep="\t", names=["seq", "name"], usecols=range(2)
        )

        # Save n calls to pd.DataFrame.__getattr__
        self.barcode_seqs = self.barcodes.seq = self.barcodes.seq.str.upper()
        # Determine the size of our design
        self.nbc = len(self.barcodes)

        # Map barcode sequences to roe indices for faster access
        self.bc_dict: dict[str, int] = {seq: i for i, seq in self.barcode_seqs.items()}
        self.barcode_seqs_l: list[str] = self.barcode_seqs.to_list()
        self.barcode_masks_l: list[tuple[int]] = [
            Demultiplexer.seq_to_mask(seq) for seq in self.barcode_seqs
        ]

        self.barcodes.loc[-1] = ("N" * self.bclen, "noBC")
        self.barcode_seqs = self.barcodes.seq
        self.barcode_names = self.barcodes.name

        # Determine the number of output files
        self.nfiles = self.split or self.nbc

        # Prefabricate python slices for quick indexing
        self.bcslice = slice(self.bclen)
        self.umislice = slice(self.umistart - 1, self.umistart + self.umilen - 1)

        # QC statistics tracker
        self.stats = ReadCounts(
            0,
            {name: 0 for name in self.barcode_names},
            [0 for _ in range(self.nfiles + 1)],
            [
                {name: 0 for name in self.barcode_names.drop(-1)}
                for _ in range(self.mismatch + 1)
            ],
            {
                reason.name: {name: 0 for name in self.barcode_names}
                for reason in QCFailReason
            },
            collections.Counter(),
        )

        # Validate that the barcodes are different enough for the input parameters.
        logger = logging.getLogger("BCcheck")
        if self.barcode_seqs.duplicated().any():
            logger.error(
                "Barcode sequences are not all unique. Cowardly refusing to proceed."
            )
            sys.exit(1)
        if self.barcode_names.duplicated().any():
            logger.error(
                "Barcode names are not all unique. Cowardly refusing to proceed."
            )
            sys.exit(1)
        if check and self.mismatch != 0:
            self.check()

        self.exit_stack = contextlib.ExitStack()

        # Rudimentary duplicate pre-screening
        self.seen: set[ReadPair] = set()

    @property
    def adapter(self):
        return self._adapter and BackAdapter(self._adapter)

    @adapter.setter
    def adapter(self, adapter):
        self._adapter = adapter

    def __enter__(self):
        return self

    def __exit__(self, *e):
        self.exit_stack.__exit__(*e)

    def init_filters(self):
        # Determine which QC filters to enable or not
        enable_homogeny_filter = min(self.max_AT, self.max_GC) < 1
        enable_homopolymer_filter = self.min_homop_len >= 5
        enable_quality_trimming = self.min_qual > 0
        enable_umi_quality_filter = (
            enable_quality_trimming and self.umilen >= self.min_qual_bases
        )
        enable_adapter_filter = self.adapter is not None and len(self.adapter) >= 5
        self.enabled_checks: list[Callable[[ReadPair], QCFailReason]] = [
            self.should_discard_length
        ]
        if enable_umi_quality_filter:
            self.enabled_checks.append(self.should_discard_umi_quality)
        if enable_homogeny_filter:
            self.enabled_checks.append(self.should_discard_homogeny)
        if enable_homopolymer_filter:
            self.enabled_checks.append(self.should_discard_homopolymer)
        if enable_adapter_filter:
            self.enabled_checks.append(self.should_discard_adapter)
        if enable_quality_trimming:
            self.enabled_checks.append(self.should_discard_trimming)
        # self.enabled_checks.append(self.should_discard_duplicate)

        # Precompute the homopolymer adapters
        self.homopolymers: frozenset[str] = frozenset(
            (c * self.min_homop_len for c in "ACGTN")
            if enable_homopolymer_filter
            else ()
        )

    def check(self):
        logger = logging.getLogger("BCcheck")
        self.mismatch *= 2
        too_close = False
        for i, barcode in enumerate(self.barcode_seqs_l[:-1]):
            idx, mismatch = self.match_barcode(barcode, search_offset=i + 1)
            if idx != -1:
                logger.error(
                    "Barcodes %d and %d are too close! (%d substitution distance)",
                    i,
                    idx,
                    mismatch,
                )
                too_close = True
        if too_close:
            sys.exit(1)
        self.mismatch //= 2

    @staticmethod
    def seq_to_mask(seq: str):
        return tuple(Demultiplexer.IUPAC_ARR[ord(c)] for c in seq)

    @cache
    def get_homopolymer_regex(self, adaptor: str):
        """
        Gets (and caches) the regex.Pattern corresponding to the adaptor with substitutions allowed
        :param adaptor: Adaptor sequence to (closely) match
        :return: A precompiled regex.Pattern representing the adaptor and up to self.mismatch substitutions.
        """
        # regex aped from STPipeline
        return regex.compile(rf"(?:{adaptor}){{s<={self.max_homop_mmatch}}}")

    def find_homopolymer(self, rp: ReadPair, adaptor: str) -> int:
        """
        Implements STPipeline's homopolymer search, allowing for mismatches
        :param rp: The ReadPair to search
        :param adaptor: The adapter sequence to query
        :return: The integer index of the adaptor in the sequence, or -1 if not found
        """
        cdna_seq = rp.cdna_seq
        if len(cdna_seq) < len(adaptor):  # Read is too short to trim
            return -1
        elif self.max_homop_mmatch == 0:  # Find exact matches only
            return cdna_seq.find(adaptor)
        elif m := self.get_homopolymer_regex(adaptor).search(
            cdna_seq
        ):  # TODO: this is slow
            return m.start() + (adaptor[0] != m[0][0] and m[0].find(adaptor[0]))
        else:  # Not found
            return -1

    def iter_read_pairs(
        self, r1fqname: str, r2fqname: str, *, phred=33
    ) -> Iterator[ReadPair]:
        with pysam.FastxFile(r1fqname, persist=False) as r1fq, pysam.FastxFile(
            r2fqname, persist=False
        ) as r2fq:
            for r1, r2 in zip(r1fq, r2fq):
                # while r1.name != r2.name:
                #     r1 = next(r1fq)
                head_r2 = r2.name
                seq_r1 = r1.sequence.upper()
                seq_r2 = r2.sequence.upper()
                qual_r1 = r1.quality
                r1qa = r1.get_quality_array(phred)
                r2qa = r2.get_quality_array(phred)
                if len(seq_r1) < self.bclen + self.umilen:
                    raise IndexError(
                        f"R1 length too short (expected {self.bclen + self.umilen}, got {len(seq_r1)})"
                    )
                yield ReadPair(
                    head_r2.split()[0],
                    seq_r2,
                    r2qa,
                    seq_r1[self.bcslice],
                    qual_r1[self.bcslice],
                    seq_r1[self.umislice],
                    qual_r1[self.umislice],
                    r1qa[self.bcslice],
                    r1qa[self.umislice],
                )

    def match_sequences(self, seqA: tuple[int], seqB: tuple[int]):
        """
        Helper function for match_barcode
        :param seqA: Integer-encoded sequence A
        :param seqB: Integer-encoded sequence B
        :return: Number of mismatches, capping out at 1 + the maximum
        """
        return sum(not a & b for a, b in zip(seqA, seqB))

    @lru_cache(1024)
    def match_barcode(self, seq: str, *, search_offset=0) -> tuple[int, int]:
        """
        Looks up the barcode sequence in the whitelist. If
        an exact match is not found, finds the closest match
        within the maximum number of mismatches.
        If multiple barcode sequences would match, the first
        in the whitelist is selected. Researchers should
        be careful to select barcode sequences and a mismatch
        tolerance that would reasonably avoid such ambiguities.

        :param seq: The barcode sequence
        :param search_offset: kw-only first index to consider, internal use only
        :return:
            A 2-tuple containing the 0-based index of the barcode sequence and the number of mismatches.
            If not found, both are -1.
        """
        candidates = self.barcode_seqs_l
        if seq.count("N") > self.max_bc_ns:
            return -1, -1
        if search_offset == 0:
            idx = self.bc_dict.get(seq, -1)
        else:
            try:
                idx = candidates.index(seq, search_offset)
            except ValueError:
                idx = -1
        if idx >= 0:
            return idx, 0
        if not self.mismatch:
            return -1, -1
        mask = Demultiplexer.seq_to_mask(seq)

        result = min(
            enumerate(
                (
                    self.match_sequences(mask, candidate)
                    for candidate in self.barcode_masks_l[search_offset:]
                ),
                search_offset,
            ),
            key=lambda t: t[1],
        )
        if result[1] > self.mismatch:
            return -1, -1
        return result

    def demultiplex_read(self, rp: ReadPair) -> tuple[int, int, int]:
        """
        Matches a read to a barcode in the whitelist. As a side effect, updates the internal demultiplexing statistics.
        :param rp: A ReadPair object i.e. from iter_read_pairs
        :return: A 3-tuple containing the barcode index, number of substitutions, and file index.
        """
        idx, mismatch = self.match_barcode(rp.bcseq)
        fileidx = idx * self.split // self.nbc if idx >= 0 and self.split else idx
        return idx, mismatch, fileidx

    def quality_trim_index(self, rp: ReadPair):
        """
        Helper function to determine the 3' base position to trim from. Adapted from STPipeline who aped it from
        CutAdapt. The license for this snippet is replicated below.

        https://github.com/marcelm/cutadapt/
        Copyright (c) 2010-2022 Marcel Martin <marcel.martin@scilifelab.se>

        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in
        all copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
        THE SOFTWARE.

        :param rp: The ReadPair to trim
        :return: The index position at the 3' end beyond which bases and qualities ought to be trimmed.
        """
        s = 0
        max_qual = 0
        max_i = len(rp.cdna_qual)
        for i in reversed(range(max_i)):
            q = rp.cdna_qual[i]
            if rp.cdna_seq[i] == "G":
                q = self.min_qual - 1
            s += self.min_qual - q
            if s < 0:
                break
            if s > max_qual:
                max_qual = s
                max_i = i
        return max_i

    def should_discard_umi_quality(self, rp: ReadPair) -> QCFailReason:
        """
        Helper function not to be called directly.

        Tests if not enough UMI bases pass the quality threshold.
        """
        return (
            QCFailReason.FAIL_UMI_QUALITY
            if sum(x < self.min_qual for x in rp.umiqual_a) > self.min_qual_bases
            else QCFailReason.PASS
        )

    def should_discard_length(self, rp: ReadPair) -> QCFailReason:
        """
        Helper function not to be called directly.

        Tests if the cDNA read is long enough to proceed.
        """
        return (
            QCFailReason.FAIL_LENGTH
            if len(rp.cdna_seq) < self.min_length
            else QCFailReason.PASS
        )

    def should_discard_homogeny(self, rp: ReadPair) -> QCFailReason:
        """
        Helper function not to be called directly.

        Tests if the AT and GC content fraction are small enough to proceed.
        """
        cdna_seq = rp.cdna_seq
        clen = len(cdna_seq)
        n_AT = n_GC = 0
        for c in cdna_seq:
            if c in "AT":
                n_AT += 1
            elif c in "CG":
                n_GC += 1
        if n_AT >= clen * self.max_AT:
            return QCFailReason.FAIL_AT_CONTENT
        elif n_GC >= clen * self.max_GC:
            return QCFailReason.FAIL_GC_CONTENT
        else:
            return QCFailReason.PASS

    def should_discard_homopolymer(self, rp: ReadPair) -> QCFailReason:
        """
        Helper function not to be called directly.

        Trims homopolymers and tests that the result is still long enough to proceed.
        """
        for homop in self.homopolymers:
            if (homo_idx := self.find_homopolymer(rp, homop)) != -1:
                clen = homo_idx
                if clen < self.min_length:
                    return QCFailReason.FAIL_HOMOPOLYMER
                rp.trim(clen)
        return QCFailReason.PASS

    def should_discard_adapter(self, rp: ReadPair) -> QCFailReason:
        """
        Helper function not to be called directly.

        Tests if trimming sequencing adapters results in a long-enough read.
        """
        match = self.adapter.match_to(rp.cdna_seq)
        if match is None:
            return QCFailReason.PASS
        if match.rstart < self.min_length:
            return QCFailReason.FAIL_ADAPTER_TRIMMING
        # rp = match.trimmed(rp)
        rp.trim(match.rstart)
        return QCFailReason.PASS

    def should_discard_trimming(self, rp: ReadPair) -> QCFailReason:
        """
        Helper function not to be called directly.

        Tests if trimming low-quality 3' bases results in a long-enough read.
        """
        trim_idx = self.quality_trim_index(rp)
        if trim_idx < self.min_length:
            return QCFailReason.FAIL_QUALITY_TRIMMING
        rp.trim(trim_idx)
        return QCFailReason.PASS

    def should_discard_duplicate(self, rp: ReadPair) -> QCFailReason:
        """
        Helper function not to be called directly

        Tests if a read with the exact same sequence, barcode, and UMI has already been observed
        """
        ret = QCFailReason.FAIL_DUPLICATE if rp in self.seen else QCFailReason.PASS
        self.seen.add(rp)
        return ret

    def quality_check(self, rp: ReadPair) -> QCFailReason:
        """
        Performs QC checks on a ReadPair, following the specifications implemented in stpipeline.
        Also trims homopolymers and low-quality 3' bases.
        :param rp: The ReadPair to QC and trim
        :return: An enum representing the reason why the read pair failed QC (or QCFailReason.PASS if it passed QC).
        """
        for method in self.enabled_checks:
            if (reason := method(rp)) is not QCFailReason.PASS:
                return reason
        return QCFailReason.PASS

    @cache
    def get_bc_seq_name(self, idx: int) -> tuple[str, str]:
        """
        Fetches the sequence and name of a barcode given its index
        :param idx: Index of the barcode
        :return: 2-tuple of str containing sequence and name
        """
        return self.barcode_seqs[idx], self.barcode_names[idx]

    def record_stats(
        self, rp: ReadPair, qcpass: QCFailReason, name: str, mismatch: int, fileidx: int
    ):
        """
        Record demultiplex stats
        :param rp: The ReadPair of interest
        :param qcpass: The reason for QC fail (or QCFailReason.PASS if passed QC)
        :param name: Name of the well
        :param mismatch: Number of substitutions between observed and actual barcode
        :param fileidx: Index of the downstream task (STAR, featureCounts, Deduplicate)
        """
        self.stats.read_count += 1
        self.stats.fail_counts[qcpass.name][name] += 1
        self.stats.written_counts[fileidx] += 1
        self.stats.well_counts[name] += 1
        if fileidx >= 0:
            self.stats.mismatch_counts[mismatch][name] += 1
        else:
            self.stats.unmatched_counts[rp.bcseq] += 1

    @cached_property
    def sam_header(self):
        # idx * self.split // self.nbc if idx >= 0 and self.split else idx
        return {
            "HD": {"VN": "1.6", "SO": "unsorted"},
            "PG": [
                {
                    "ID": os.path.splitext(os.path.basename(__file__))[0],
                    "PN": os.path.splitext(os.path.basename(__file__))[0],
                    "VN": __version__,
                    "CL": shlex.join(sys.argv),
                }
            ],
            "RG": [
                {
                    "ID": "noBC" if i < 0 else str(i).zfill(4),
                    "BC": (
                        "-".join(
                            self.barcode_seqs[
                                i
                                * self.nbc
                                // self.split : (i + 1)
                                * self.nbc
                                // self.split
                            ]
                        )
                        if i >= 0 and self.split
                        else self.barcodes.at[i, "seq"]
                    ),
                    "PL": "ILLUMINA",
                }
                for i in range(-1, self.nfiles)
            ],
        }

    @cache
    def get_writer(
        self, prefix, bcidx: int = None, qcfail: QCFailReason = QCFailReason.PASS
    ):
        if qcfail is QCFailReason.PASS:
            suffix = "_qcFail"
        elif bcidx is None:
            suffix = ""
        elif bcidx == -1:
            suffix = "_noBC"
        else:
            suffix = f"_{bcidx:04d}"
        return self.exit_stack.enter_context(
            SamWriter(f"{prefix}{suffix}.bam", self.sam_header, cpus=self.cpus)
        )

    def demultiplex_experiment(self, r1fq: str, r2fq: str, prefix):
        """
        Ingests a pair of FASTQ files from a DRUG-seq experiment, performs QC filtering, and writes each to an output
            SAM file.
        :param r1fq: File containing R1 reads (BC+UMI)
        :param r2fq: File containing R2 reads (cDNA)
        :param prefix: Output filename prefix
        """
        logger = logging.getLogger(prefix)
        try:
            logger.info("Begin")
            for i, rp in enumerate(self.iter_read_pairs(r1fq, r2fq), 1):
                idx, mismatch, fileidx = self.demultiplex_read(rp)
                qcpass = self.quality_check(rp)
                bcseq, name = self.get_bc_seq_name(idx)
                self.record_stats(rp, qcpass, name, mismatch, fileidx)
                if qcpass is QCFailReason.PASS or self.write_qcfail:
                    writer = self.get_writer(
                        prefix, fileidx if self.BAM_WRITE_MULTIPLE else None, qcpass
                    )
                    writer.write(rp, qcpass, bcseq, mismatch, fileidx)
                if i % 1000000 == 0:
                    logger.info("Processed %d reads...", i)
            logger.info("End, processed %d reads", i)
        except Exception:
            logger.critical("Aborting", exc_info=True)
            sys.exit(1)


class CLI(argparse.Namespace):
    R1: str
    R2: str
    barcodes: typing.TextIO
    prefix: str
    bclen: int = 16
    umistart: int = 17
    umilen: int = 10
    split: int = 0
    mismatch: int = 2
    min_length: int = 20
    min_qual: int = 20
    min_qual_bases: int = 6
    max_AT: float = 0.9
    max_GC: float = 0.9
    min_homop_len: int = 10
    max_homop_mmatch: int = 0
    max_bc_ns: int = 2
    check: bool = False
    cpus: int = 1
    adapter: str = None
    write_qcfail: bool = True

    _parser = argparse.ArgumentParser()
    _parser.add_argument("R1", help="Path to read-1 fastq file (barcode+UMI)")
    _parser.add_argument("R2", help="Path to read-2 fastq file (cDNA)")
    _parser.add_argument(
        "barcodes",
        type=argparse.FileType(),
        help="Path to TSV file defining well names and barcode sequences",
    )
    _parser.add_argument("prefix", help="Output file prefix")
    _parser.add_argument(
        "--bclen",
        type=int,
        default=16,
        help="Number of bases in the well barcode (default: %(default)d)",
    )
    _parser.add_argument(
        "--umistart",
        type=int,
        default=17,
        help="Start position of the UMI in the R1 read (default: %(default)d)",
    )
    _parser.add_argument(
        "--umilen",
        type=int,
        default=10,
        help="Number of bases in the UMI sequence (default: %(default)d)",
    )
    _parser.add_argument(
        "--split",
        type=int,
        default=0,
        help="Maximum number of BAM files to write (default: one per well)",
    )
    _parser.add_argument(
        "--mismatch",
        type=int,
        default=2,
        help="Maximum number of mismatches (substitutions) between observed and defined barcode sequences "
        "(default %(default)d)",
    )
    _parser.add_argument(
        "--min-length",
        type=int,
        default=20,
        help="Minimum cDNA length to pass QC, before or after trimming (default: %(default)d)",
    )
    _parser.add_argument(
        "--min-qual",
        type=int,
        default=20,
        help="Phred quality score threshold for UMI filtering and cDNA trimming (default: %(default)d)",
    )
    _parser.add_argument(
        "--min-qual-bases",
        type=int,
        default=6,
        help="Minimum number of UMI bases with phred quality score passing --min-qual (default: %(default)d)",
    )
    _parser.add_argument(
        "--max-AT",
        type=float,
        default=0.9,
        help="Maximum fraction of cDNA bases allowed to be A or T per read (default: %(default)g)",
    )
    _parser.add_argument(
        "--max-GC",
        type=float,
        default=0.9,
        help="Maximum fraction of cDNA bases allowed to be G or C per read (default: %(default)g)",
    )
    _parser.add_argument(
        "--min-homop-len",
        type=int,
        default=10,
        help="Trim the cDNA when I encounter a run of at least this many of the same nucleotide (homopolymer) "
        "(default: %(default)d)",
    )
    _parser.add_argument(
        "--max-homop-mmatch",
        type=int,
        default=0,
        help="Allow at most this many substitutions when detecting homopolymers (default: %(default)d)",
    )
    _parser.add_argument(
        "--max-bc-ns",
        type=int,
        default=2,
        help="Maximum number of N base calls allowed in the barcode sequence (default: %(default)d)",
    )
    _parser.add_argument(
        "--check-bcs",
        dest="check",
        action="store_true",
        default=False,
        help="If set, will make sure that the barcodes are unique, and the list of barcodes is compatible with the "
        "--mismatch setting.",
    )
    _parser.add_argument("--cpus", type=int, default=1, help="Number of BAM IO threads")
    _parser.add_argument(
        "--adapter", type=str, help="Sequencing adapter to trim from the 3' end of R2"
    )
    _parser.add_argument(
        "--discard-failed",
        dest="write_qcfail",
        action="store_false",
        default=True,
        help="If set, suppresses outputting QCFail reads. Useful for speed but may hinder scientific diagnostics.",
    )

    def __init__(self, args=None):
        self.__class__._parser.parse_args(args, self)

    def main(self):
        logging.basicConfig(
            level=logging.INFO,
            format="[%(asctime)s] %(levelname)s:%(name)s:%(message)s",
        )
        prefix = self.prefix
        demultiplexer = Demultiplexer(
            self.barcodes,
            bclen=self.bclen,
            umistart=self.umistart,
            umilen=self.umilen,
            split=self.split,
            mismatch=self.mismatch,
            min_length=self.min_length,
            min_qual=self.min_qual,
            min_qual_bases=self.min_qual_bases,
            max_AT=self.max_AT,
            max_GC=self.max_GC,
            min_homop_len=self.min_homop_len,
            max_homop_mmatch=self.max_homop_mmatch,
            check=self.check,
            max_bc_ns=self.max_bc_ns,
            cpus=self.cpus,
            adapter=self.adapter,
        )
        demultiplexer.demultiplex_experiment(self.R1, self.R2, prefix)
        demultiplexer.stats.unmatched_counts = dict(
            demultiplexer.stats.unmatched_counts.most_common()
        )
        with open(f"{prefix}.demux.json", "w") as ofp:
            json.dump(dataclasses.asdict(demultiplexer.stats), ofp, indent=2)


if __name__ == "__main__":
    CLI().main()
