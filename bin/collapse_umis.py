#!/usr/bin/env python

"""
Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
This file is part of DRAGoN.

This source code is licensed under the MIT License found in the
LICENSE file in the root directory of this source tree.
"""

import argparse
import collections
import contextlib
import enum
import functools
import itertools
import logging
import operator
import os
import shlex
import sys
import threading
from collections.abc import Container, Generator, Iterable, Iterator
from typing import Any, Literal, NamedTuple, Optional, Type

# Python api invokes isinstance which is slow. Bypass it and use the C class.
try:  # since biopython-1.82
    from Bio.Align._pairwisealigner import PairwiseAligner
except ImportError:
    from Bio.Align._aligners import PairwiseAligner

import numpy as np
import pandas as pd
import pysam
import scipy.io as sio
import scipy.sparse as sp
import xopen
from version import __version__

Strand = Literal["+", "-"]
ChromKey = tuple[int, Strand]


class ReadCountPrinter(threading.Thread):
    def __init__(self):
        super().__init__(target=self.run_printer)
        self._counter = 0
        self._lock = threading.Lock()
        self._logger = logging.getLogger("ReadCount")
        self._enabled = False
        self._condition = threading.Condition(lock=self._lock)

    def increment(self):
        with self._lock:
            self._counter += 1

    def get(self):
        with self._lock:
            return self._counter

    def reset(self):
        with self._lock:
            self._counter = 0

    def run_printer(self):
        with self._condition:
            while not self._condition.wait_for(lambda: not self._enabled, 30):
                self._logger.info("Processed %d reads", self._counter)

    def stop(self):
        self._enabled = False
        with self._condition:
            self._condition.notify_all()

    def __enter__(self):
        self._enabled = True
        self.start()
        return self

    def __exit__(self, exc_type, exc, exc_tb):
        self.stop()
        self.join()


class GFFRecord(NamedTuple):
    seqname: str
    source: str
    feature: str
    start: int
    end: int
    name: str
    strand: str
    score: float
    attrs: dict[str, str]

    @classmethod
    def from_file(
        cls: "Type[GFFRecord]",
        filename: str,
        *,
        feature: Optional[str] = None,
        columns: Optional[Container[str]] = None,
    ) -> "Iterable[GFFRecord]":
        """
        Reads a GTF file and yields all records of a given feature type as
        namedtuples of type GFFRecord, which has the following elements:
            seqname: Name of the chromosome the feature is on
            source: Annotation source for the feature (eg. HAVANA)
            feature: Type of feature (gene, transcript, exon, etc.)
            start: Start coordinate of the feature, 1-based
            end: End coordinate of the feature, 1-based
            name: Name of the feature, usually blank
            strand: Strand the feature is found on, usually + or -
            score: Score of the feature, not used in transcriptome annotations
            attrs: dict of gtf attribute keys to values

        :param filename: Path to the GTF file to read
        :param feature: Select only records matching this feature type. If None (the default), select all records.
        :param columns: Select only these attributes for the .attrs field. If None (the default), select all attributes.
        :return: Iterator yielding GFFRecord
        """

        with xopen.xopen(filename) as fp:
            for line in fp:
                line = line.rstrip()
                # Skip header
                if line.startswith("#"):
                    continue
                fields: list[Any] = line.split("\t")
                # Skip unwanted feature type
                if feature is not None and fields[2] != feature:
                    continue
                # Build row object
                fields[3] = int(fields[3])
                fields[4] = int(fields[4])
                fields[5] = "" if fields[5] == "." else fields[5]
                fields[6] = "" if fields[6] == "." else fields[6]
                fields[7] = float("nan" if fields[7] == "." else fields[7])
                fields[8] = {
                    key: value.strip('"')
                    for key, value in (
                        field.split(" ", 1) for field in fields[8].split("; ")
                    )
                    if columns is None or key in columns
                }
                yield cls._make(fields)


class GenomicPos(NamedTuple):
    r_id: str
    r_pos: int
    r_end: int
    r_strand: Strand


class FilterKey(NamedTuple):
    """
    namedtuple whose attributes identify a read for the sake of duplicate filtering
        cb: Barcode sequence
        q_name: Original name of the read
        coords: Dict mapping the coordinates of all alignments, for quick overlap checking
        mapq: dict quality, either 255 for unique reads or int(-10*log10(1-1/Nmap)) for multimappers,
            as specified by STAR
        umi: UMI sequence
        gene_idxs: frozenset of indices into the set of annotated genes, to which this read was assigned
        mx: tuple of GenomicPos representing all alignments in order of appearance
        xt: tuple of frozenset of gene names representing all genes assigned to each alignment in order of appearance
    """

    cb: str
    q_name: str
    coords: dict[ChromKey, tuple[int]]
    mapq: int
    umi: str
    gene_idxs: frozenset[int]
    mx: tuple[GenomicPos]
    xt: tuple[frozenset[str]]


class DuplicateFiltrator:
    def __init__(
        self,
        header: pysam.AlignmentHeader,
        umilen: int,
        max_distance: int,
        max_lookback: int,
        aligner: PairwiseAligner,
    ):
        """
        Creates a DuplicateFiltrator. One of these will be spawned for each set of genes.
        :param umilen: Length of the UMI
        :param min_score: Minimum alignment score between UMI sequences
        :param max_lookback: Maximum genomic distance between two UMIs
        :param aligner: Reference to a PairwiseAligner instance
        """
        self.memory: dict[ChromKey, dict[int, dict[str, str]]] = {
            (i, s): {} for i in range(header.nreferences) for s in "+-"
        }
        self.umilen = umilen
        self.max_distance = max_distance
        self.min_score = umilen - max_distance
        self.max_lookback = max_lookback
        self.score = functools.lru_cache(1000000)(aligner.score)

    def __call__(self, umi: FilterKey) -> Optional[str]:
        """
        Determines if the read is a duplicate based on what's already been seen.
        :param umi: A FilterKey representing the UMI to check
        :return: If the input UMI is a duplicate, returns the name of the read it is a duplicate of.
        """
        for chrom, posns in umi.coords.items():
            if not (chrom_d := self.memory.get(chrom)):
                continue
            for pos in posns:
                for i in range(pos - self.max_lookback, pos + self.max_lookback + 1):
                    if (pos_l := chrom_d.get(i)) is None:
                        continue
                    if second := pos_l.get(umi.umi):
                        return second
                    if self.max_distance == 0:
                        continue
                    for umi2, second in pos_l.items():
                        if self.score(umi.umi, umi2, "+") >= self.min_score:
                            return second
        # New UMI, stick it in memory
        for chrom, posns in umi.coords.items():
            chrom_d = self.memory.setdefault(chrom, {})
            for pos in posns:
                chrom_d.setdefault(pos, {})[umi.umi] = umi.q_name

    def reset(self):
        """
        Reset the underlying deque
        """
        self.memory.clear()


class MultiMapStrategy(enum.Enum):
    """
    String enum class encoding the strategy for distributiing multimapped UMIs
    """

    UNIQUE = "Unique"  # Discard multimapped UMIs
    UNIFORM = "Uniform"  # Distribute multimapped UMIs evenly across all loci
    PROPUNIQUE = "PropUnique"  # Distribute UMIs proportionately to the unique counts
    EM = "EM"  # Expectation maximization
    RESCUE = "Rescue"  # The first iteration of EM


class Stat(enum.Enum):
    """
    String enum class encoding the columns of the stats report
    """

    TOTAL = "total"
    DISCARD_QCFAIL = "discard_qcFail"
    DISCARD_UNMAPPED = "discard_unmapped"
    DISCARD_NOFEATURE = "discard_noFeature"
    PASSED_QC = "passed_qc"
    N_MULTIMAPPED = "n_multimapped"
    UMI_COUNTS = "umi_counts"
    UNIQUE_COUNTS = "unique_counts"
    UNIQUE_UMIS = "unique_umis"


class Deduplicator:
    def __init__(
        self,
        bamfilename: str,
        threads: int,
        barcodes_in: str,
        gtf_in: str,
        umilen: int = 10,
        max_distance: int = 2,
        max_lookback: int = 10,
        mmap_strategies: frozenset[MultiMapStrategy] = frozenset(
            (MultiMapStrategy.UNIQUE,)
        ),
        em_iter_max: int = 100,
        em_max_diff: float = 0.01,
        em_small_thresh: float = 0.01,
        keep_bam: bool = False,
    ):
        """
        Wrapper class for removing PCR duplicates (using the UMI as a backbone).
        :param bamfilename: Path to BAM file containing reads demultiplexed by demultiplex.py, mapped by STAR,
                        and assigned to genes by featureCounts vis-a-vis assign_featurecounts.sh.
        :param threads: Number of parallel threads for IO
        :param barcodes_in: Path to the barcodes TSV
        :param gtf_in: Path to the GTF file
        :param umilen: Length of the UMI sequence
        :param max_distance: Maximum edit distance between two UMIs to consider them the same (default: 2)
        :param max_lookback: Maximum genomic distance to test for matching UMIs (default: 10)
        :param mmap_strategies: A MultiMapStrategy (default: MultiMapStrategy.UNIQUE)
        :param em_iter_max: If mmap_strategies contains MultiMapStrategy.EM, the maximum number of iterations of EM
                        (default: 100)
        :param em_max_diff: If mmap_strategies contains MultiMapStrategy.EM, the tolerance for convergence, measured by
                        the largest absolute change in counts for any read in the current well (default: 0.01)
        :param em_small_thresh: If mmap_strategies contains MultiMapStrategy.EM, clamp counts smaller than this to 0
                        (default: 0.01)
        :param keep_bool: If True, outputs a BAM file
        """
        self.logger = logging.getLogger("Deduplicate")
        self.logger.info("Begin deduplicate workflow")
        self.logger.info("Reading barcodes %s", barcodes_in)
        self.barcodes = pd.read_csv(
            barcodes_in, sep="\t", names=["seq", "well"], usecols=range(2)
        )
        self.logger.info("Reading features %s", gtf_in)
        self.features = pd.DataFrame(
            [
                [
                    row.attrs.get("gene_id", ""),
                    row.attrs.get("gene_name", ""),
                    "Gene Expression",
                ]
                for row in GFFRecord.from_file(
                    gtf_in, feature="gene", columns=("gene_id", "gene_name")
                )
            ],
            columns=["Geneid", "GeneName", "FeatureType"],
        )
        self._context_stack = contextlib.ExitStack()
        self.bamfilename = bamfilename
        self.iothreads = threads
        self.bamfile: Optional[pysam.AlignmentFile] = None
        self.outbam: Optional[pysam.AlignmentFile] = None
        self.umilen = umilen
        self.max_distance = max_distance
        self.max_lookback = max_lookback
        self.keep_bam = keep_bam

        # In the default execution, this will only run on one barcode (well).
        # However, the user can group barcodes using a nextflow parameter.
        # Also, the counts matrix specification is to indicate the index of the barcode
        self.matrices = {
            x: np.zeros(
                (len(self.features), len(self.barcodes)),
                dtype=int if x is MultiMapStrategy.UNIQUE else float,
            )
            for x in mmap_strategies
        }
        self.counts_per_read: dict[MultiMapStrategy, sp.lil_array] = {}
        self.umi_names: dict[str, int] = {}

        # Hashmaps for quick name-to-index lookup
        self.barcode_seqs = {seq: i for i, seq in enumerate(self.barcodes.seq)}
        self.barcode_names = {
            seq: name for seq, name in zip(self.barcodes.seq, self.barcodes.well)
        }
        self.feature_names = {name: i for i, name in enumerate(self.features.Geneid)}
        self.geneids_l = list(self.features.Geneid)

        # Precompute some frequently-used attributes and preallocate containers
        self.stats = {
            stat: {name: 0 for name in self.barcodes.well.to_list() + ["noBC"]}
            for stat in Stat
        }
        self.aligner = PairwiseAligner()
        self.aligner.gap_score = -0.5  # Explicitly penalize indels
        self.min_score = self.umilen - self.max_distance
        self.multimappers: collections.Counter[frozenset[int]] = collections.Counter()
        self.mmap_strategies = sorted(
            mmap_strategies - {MultiMapStrategy.UNIQUE},
            key=list(MultiMapStrategy).index,
        )
        self.is_counting_multimappers = bool(self.mmap_strategies)
        self.em_iter_max = em_iter_max
        self.em_max_diff = em_max_diff
        self.em_small_thresh = em_small_thresh

        self.printer = ReadCountPrinter()

    def __enter__(self):
        self.bamfile = self._context_stack.enter_context(
            pysam.AlignmentFile(self.bamfilename, threads=self.iothreads)
        )
        header = self.bamfile.header.to_dict()
        pg = {
            "ID": os.path.splitext(os.path.basename(__file__))[0],
            "PN": os.path.splitext(os.path.basename(__file__))[0],
            "VN": __version__,
            "CL": shlex.join(sys.argv),
        }
        if "PG" in header:
            header["PG"].append(pg)
        else:
            header["PG"] = [pg]
        stem, ext = os.path.splitext(self.bamfilename)
        if self.keep_bam:
            self.outbam = self._context_stack.enter_context(
                pysam.AlignmentFile(
                    f"{stem}_MarkDups{ext}",
                    "w" if ext == ".sam" else "wbu",
                    header=header,
                    threads=self.iothreads,
                )
            )
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._context_stack.__exit__(exc_type, exc_val, exc_tb)
        self.outbam = None
        self.bamfile = None

    def filter_reads(self) -> Generator[tuple[FilterKey, bool, pysam.AlignedSegment]]:
        """
        Reads the BAM file and discards any reads failing certain criteria.
        Assumes the BAM file is sorted by query name.
        :return: A generator yielding 3-tuples of FilterKey, boolean discard flag, and the original aligned segment
        """

        # Previously collapse_multimappers.py
        group_key = operator.attrgetter("query_name")
        sort_key = operator.attrgetter("is_secondary")

        group: Iterable[pysam.AlignedSegment]
        for name, group in itertools.groupby(self.bamfile, group_key):
            primary, *secondaries = group = sorted(group, key=sort_key)
            if primary.is_secondary:
                print("Sorting by name did not work as expected", file=sys.stderr)
                sys.exit(1)
            cb = primary.get_tag("CB")
            mx = tuple(
                GenomicPos(
                    read.reference_name,
                    read.reference_start,
                    read.reference_end,
                    "+-"[read.is_reverse],
                )
                for read in group
            )
            sx: collections.defaultdict[ChromKey, tuple[int]] = collections.defaultdict(
                tuple
            )
            for pos in mx:
                sx[(pos.r_id, pos.r_strand)] += (pos.r_pos,)
            xt = [read.get_tag("XT") if read.has_tag("XT") else "-" for read in group]
            gene_idxs = frozenset(
                self.feature_names.get(g) for f in xt if f != "-" for g in f.split(",")
            )
            key_xt = tuple(frozenset(xti.split(",")) for xti in xt)
            if secondaries:
                xs = [read.get_tag("XS") for read in group]
                primary.set_tag(
                    "MX",
                    ";".join(
                        f"{r.r_id}{r.r_strand}{r.r_pos}-{r.r_end}" for r in mx[1:]
                    ),
                )
                if gene_idxs:
                    primary.set_tag("XN", len(gene_idxs))
                primary.set_tag("XS", ";".join(xs))
                primary.set_tag("XT", ";".join(xt))
            discard = False

            # Record stats
            bcname = self.barcode_names.get(cb, "noBC")
            self.stats[Stat.TOTAL][bcname] += 1  # Count the read as it comes in
            if bcname == "noBC":
                discard = True
            if (
                primary.is_qcfail or primary.get_tag("QC") != "PASS"
            ):  # read excluded because QCFail
                self.stats[Stat.DISCARD_QCFAIL][bcname] += 1
                discard = True
            elif primary.is_unmapped:  # reads mapping to no genomic site
                self.stats[Stat.DISCARD_UNMAPPED][bcname] += 1
                discard = True
            elif not gene_idxs:  # reads not overlapping any feature
                self.stats[Stat.DISCARD_NOFEATURE][bcname] += 1
                discard = True
            else:
                self.stats[Stat.PASSED_QC][bcname] += 1
                if secondaries:
                    self.stats[Stat.N_MULTIMAPPED][bcname] += 1

            # Yield the results
            key = FilterKey(
                cb,
                name,
                sx,
                primary.mapping_quality,
                primary.get_tag("UB"),
                gene_idxs,
                mx,
                key_xt,
            )
            yield key, discard, primary

    def _handle_multimappers_Uniform(self, rg: int, multimaps: dict[tuple[int], int]):
        """
        Should not be called directly

        Uniformly distributes multi-mapped UMIs
        """
        uniques = self.matrices[MultiMapStrategy.UNIQUE][:, rg].astype(float)
        uniform = uniques.copy()
        for key, reads in multimaps.items():
            uniform[key,] += reads / len(key)
        self.matrices[MultiMapStrategy.UNIFORM][:, rg] = uniform

    def _handle_multimappers_PropUnique(
        self, rg: int, multimaps: dict[tuple[int], int]
    ):
        """
        Should not be called directly

        Distributes multi-mapped UMIs proportional to uniquely-mapped UMIs
        """
        uniques = self.matrices[MultiMapStrategy.UNIQUE][:, rg].astype(float)
        prop_unique = uniques.copy()
        for key, reads in multimaps.items():
            of_interest = uniques[key,]
            tot = of_interest.sum()
            if tot == 0:
                # No unique mappers, distribute uniformly
                contribution = 1 / len(key)
            else:
                # Distribute proportionally
                contribution = of_interest / tot
            prop_unique[key,] += contribution * reads
        self.matrices[MultiMapStrategy.PROPUNIQUE][:, rg] = prop_unique

    def _handle_multimappers_EM(self, rg: int, multimaps: dict[tuple[int], int]):
        """
        Should not be called directly

        Uses EM to distribute multimapped UMIs to genes
        """
        uniques = self.matrices[MultiMapStrategy.UNIQUE][:, rg].astype(float)
        # Start with sum of unique and uniform
        em = uniques.copy()
        iter_max = self.em_iter_max
        max_diff = self.em_max_diff
        small_thresh = self.em_small_thresh
        for key, reads in multimaps.items():
            em[key,] += reads / len(key)
        for _ in range(iter_max):
            em_old = em
            em_old[em_old < small_thresh] = 0  # zero-out very small counts
            em = uniques.copy()
            for key, reads in multimaps.items():
                tmp = em_old[key,]
                em[key,] += tmp * reads / tmp.sum()
            if np.abs(em - em_old).max() < max_diff:
                break
        self.matrices[MultiMapStrategy.EM][:, rg] = em

    def _handle_multimappers_Rescue(self, rg: int, multimaps: dict[tuple[int], int]):
        """
        Should not be called directly

        Uses the first round of EM to distribute multimapped UMIs to genes
        """
        uniques = self.matrices[MultiMapStrategy.UNIQUE][:, rg].astype(float)
        rescue = uniques.copy()
        for key, reads in multimaps.items():
            of_interest = uniques[key,]
            tot = of_interest.sum()
            contribution = reads * (of_interest + reads / len(key)) / (tot + reads)
            rescue[key,] += contribution
        self.matrices[MultiMapStrategy.RESCUE][:, rg] = rescue

    def handle_multimappers(self, rg: int):
        """
        Wrapper function for handling multimapped UMIs.
        :param rg: Name of the well corresponding to the current barcode
        :param multimaps: Mapper from frozenset of genes to number of UMIs
        """
        logger = logging.getLogger("handle_multimappers")
        multimaps = {tuple(key): value for key, value in self.multimappers.items()}
        strats = self.mmap_strategies
        for strat in strats:
            logger.debug("strat: %s", strat.name)
            getattr(self, f"_handle_multimappers_{strat.value}")(rg, multimaps)

    @functools.cached_property
    def filtrator(self):
        """
        Generates a DuplicateFiltrator from this class instance
        :return: A DuplicateFiltrator object
        """
        return DuplicateFiltrator(
            self.bamfile.header,
            self.umilen,
            self.max_distance,
            self.max_lookback,
            self.aligner,
        )

    def process_cb_umis(
        self, cb: str, group: Iterable[tuple[FilterKey, bool, pysam.AlignedSegment]]
    ):
        """
        Intermediate iterator to propeerly handle multimapping reads
        :param cb: Barcode sequence
        :param group: Filtered reads generator
        """
        logger = logging.getLogger("process_cb_umis")

        self.filtrator.reset()
        self.multimappers.clear()
        cb_idx = self.barcode_seqs.get(cb, -1)
        cb_name = self.barcode_names.get(cb, -1)

        key: FilterKey
        discard: bool
        read: pysam.AlignedSegment
        self.umi_names.clear()
        umi_i = 0

        logger.info("Begin processing UMIs for barcode %s", cb)

        for key, discard, read in group:  # group gives all reads in the current well
            self.printer.increment()
            if not discard:  # a potentially valid UMI
                if is_unique := len(key.gene_idxs) == 1:
                    self.stats[Stat.UNIQUE_COUNTS][cb_name] += 1
                if (
                    orig := self.filtrator(key)
                ) is None:  # callable class to filter UMIs
                    read.is_duplicate = False  # clear this BAM flag if set
                    self.stats[Stat.UMI_COUNTS][cb_name] += 1
                    if is_unique:
                        # UMI is uniquely mapped, count it now
                        (gene,) = key.gene_idxs
                        self.matrices[MultiMapStrategy.UNIQUE][gene, cb_idx] += 1
                        self.stats[Stat.UNIQUE_UMIS][cb_name] += 1
                    elif self.is_counting_multimappers:
                        self.umi_names[key.q_name] = umi_i
                        umi_i += 1
                        self.multimappers[key.gene_idxs] += 1
                elif self.keep_bam:
                    read.is_duplicate = True
                    read.set_tag("DN", orig)
                logger.debug(
                    "Well %d, Read %d, UMI %d",
                    cb_idx,
                    self.stats[Stat.PASSED_QC][cb_name],
                    self.stats[Stat.UMI_COUNTS][cb_name],
                )
            elif self.keep_bam:
                read.is_qcfail = True
            if self.keep_bam:
                self.outbam.write(read)
        # We've processed all reads assigned to this well, so let's apply the multimap strategy
        if self.is_counting_multimappers:
            logger.info("Handling multimaps for barcode %s", cb)
            self.handle_multimappers(cb_idx)
        logger.info("Finished barcode %s with %d reads", cb, self.printer.get())
        self.printer.reset()

    def get_umis(self):
        """
        The main workhorse of this modue, streams reads from the BAM file, which is sorted
        by read group and coordinate. Calls PCR duplicates by accepting the
        first read mapping to a given locus.
        :return: self
        """
        # Assumes sorted by read group, see map_star.sh#L19
        # To simulate a sliding window
        # This may look like a nested for loop, but it's linear in the number of SAM reads.
        # See documentation for itertools.groupby
        with self.printer:
            for key, group in itertools.groupby(self.filter_reads(), lambda t: t[0].cb):
                self.process_cb_umis(key, group)
        return self

    def dump_matrix(self, path: str):
        """
        Exports the counts matrix in 10x format to the specified directory.
        Writes the barcode names to barcodes.tsv, gene names to features.tsv,
        and counts matrix in sparse format to {mtx_name}.mtx.
        :param path: Path to output directory
        :return: self
        """
        os.makedirs(path, exist_ok=True)
        self.barcodes.seq.to_csv(f"{path}/barcodes.tsv", index=False, header=False)
        self.features.to_csv(
            f"{path}/features.tsv", index=False, header=False, sep="\t"
        )
        for strat, matrix in self.matrices.items():
            mtx_name = (
                "matrix"
                if strat is MultiMapStrategy.UNIQUE
                else f"UniqueAndMult-{strat.value}"
            )
            sio.mmwrite(f"{path}/{mtx_name}.mtx", sp.csc_matrix(matrix))
        return self


class CLI(argparse.Namespace):
    bamfile: str
    barcodes: str
    annofile: str
    outdir: str
    cpus: int = 1
    umilen: int = 10
    max_distance: int = 2
    max_lookback: int = 10
    mmap_strategies: list[str] = ["Unique"]
    em_iter_max: int = 100
    em_max_diff: float = 0.01
    em_small_thresh: float = 0.01
    debug: bool = False
    keep_bam: bool = False

    _parser = argparse.ArgumentParser()
    _parser.add_argument("bamfile", help="Annotated BAM file")
    _parser.add_argument(
        "barcodes",
        help="TSV file mapping barcode sequences to well names and coordinates",
    )
    _parser.add_argument("annofile", help="GTF file defining quantified features")
    _parser.add_argument("outdir", help="Output directory")
    _parser.add_argument(
        "--cpus",
        type=int,
        default=1,
        help="Number of threads for BAM IO (default: %(default)d)",
    )
    _parser.add_argument(
        "--umilen",
        type=int,
        default=10,
        help="UMI sequence length (default: %(default)d)",
    )
    _parser.add_argument(
        "--max-distance",
        type=int,
        default=2,
        help="Maximum edit distance between two UMIs to consider them as duplicates "
        "(default: %(default)d)",
    )
    _parser.add_argument(
        "--max-lookback",
        type=int,
        default=10,
        help="Maximum genomic distance between two UMIs to test them as duplicates "
        "(default: %(default)d)",
    )
    _parser.add_argument(
        "--mmap-strategies",
        nargs="+",
        default=[MultiMapStrategy.UNIQUE.value],
        choices=[x.value for x in MultiMapStrategy],
        help="Strategy(ies) for handling reads assigned to multiple features "
        "(default: %(default)s)",
    )
    _parser.add_argument(
        "--em-iter-max",
        type=int,
        default=100,
        help="If --ambiguous contains EM, the maximum number of iterations of EM "
        "(default: %(default)d)",
    )
    _parser.add_argument(
        "--em-max-diff",
        type=float,
        default=0.01,
        help="If --ambiguous contains EM, the tolerance for convergence, measured by the largest "
        "absolute change in counts for any gene in the current well (default: %(default).2f)",
    )
    _parser.add_argument(
        "--em-small-thresh",
        type=float,
        default=0.01,
        help="If --ambiguous contains EM, clamp counts smaller than this to 0 "
        "(default: %(default).2f)",
    )
    _parser.add_argument(
        "--debug", action="store_true", default=False, help="Increase logging verbosity"
    )
    _parser.add_argument(
        "--keep-bam",
        action="store_true",
        default=False,
        help="If set, outputs a BAM file called {stem}_MarkDups.bam",
    )

    def __init__(self, args=None):
        self.__class__._parser.parse_args(args, self)

    def main(self):
        logging.basicConfig(
            level=logging.DEBUG if self.debug else logging.INFO,
            format="[%(asctime)s] %(levelname)s:%(name)s:%(message)s",
        )
        with Deduplicator(
            self.bamfile,
            self.cpus,
            self.barcodes,
            self.annofile,
            umilen=self.umilen,
            max_distance=self.max_distance,
            max_lookback=self.max_lookback,
            mmap_strategies=frozenset(
                MultiMapStrategy(x) for x in self.mmap_strategies
            ),
            em_iter_max=self.em_iter_max,
            em_max_diff=self.em_max_diff,
            em_small_thresh=self.em_small_thresh,
            keep_bam=self.keep_bam,
        ) as dedup:
            dedup.get_umis()
        dedup.dump_matrix(self.outdir)
        pd.DataFrame({s.name: v for s, v in dedup.stats.items()}).to_csv(
            f"{self.outdir}/dedup_stats.tsv", sep="\t"
        )


if __name__ == "__main__":
    CLI().main()
