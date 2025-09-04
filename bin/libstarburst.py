# Library of function shared across Python scripts
# Scott Norton
# Version 2025.3.18

"""
Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
This file is part of TDP-DISCOVERY-INFORMATICS-DRUG-SEQ-PIPELINE.

This source code is licensed under the MIT License found in the
LICENSE file in the root directory of this source tree.
"""

import argparse
import collections.abc
import functools
import logging
import os
import re
import typing
import uuid

try:
    from multidict import MultiDict as mdict
except ImportError:
    mdict = dict


class GtfParseError(Exception):
    pass


def wrap_exception(
    catch_exc: type[BaseException] | tuple[type[BaseException]],
    wrap_exc: type[BaseException],
    *exc_args,
    **exc_kwargs,
):
    def wrapper(func):
        @functools.wraps(func)
        def inner(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except catch_exc as e:
                raise wrap_exc(*exc_args, **exc_kwargs) from e

        return inner

    return wrapper


def logs_runtime(func):
    """
    Logs start and finish times for the wrapped process.
    Will create a logger with a unique ID for each call to the wrapped function
    Pass a logger via the `logger` kwarg to the wrapped function to use that instead
    """
    logger = logging.getLogger(f"{func.__name__}:{uuid.uuid4().int % 1_000_000_000}")

    @functools.wraps(func)
    def inner(*args, **kwargs):
        my_logger: logging.Logger = kwargs.pop("logger", logger)
        my_logger.info("Begin")
        ret = func(*args, **kwargs)
        my_logger.info("Finish")
        return ret

    return inner


gtf_row = tuple[str, str, str, int, int, float, str, int | None, mdict[str, str]]


@wrap_exception((TypeError, ValueError, argparse.ArgumentTypeError), GtfParseError)
@logs_runtime
def read_gtf(
    fname: str | bytes | os.PathLike | typing.TextIO,
) -> collections.abc.Generator[gtf_row, typing.Any, None]:
    """Returns a generator where each row represents a GTF

    Args:
        fname (str | bytes | os.PathLike | typing.TextIO): Filename or handle.
                If a path is passed, it will be opened read-only and closed after.

    Yields:
        tuple[str, str, str, int, int, float, str, int, Mapping]:
            - seqname: Chromosome where the feature is mapped
            - source: Annotation source (ENSEMBL, Havana, etc.)
            - feature: Feature type tag (gene, transcript, exon)
            - start: 0-based start coordinate
            - end: 0-based end coordinate
            - score: Feature score
            - strand: Strand (+ or -) where the feature is mapped
            - frame: For exons, the codon offset of the first position within the CDS
            - attr: A mapping from attribute key to values. If multidict is installed, this
                    is a `class:multidict.Multidict`. Otherwise, this is a `type:dict`.
    """

    def parse_attrs(val):
        return mdict(
            (m[1], m[2].strip('"'))
            for m in re.finditer(r'(\w+) ("[^"]+"|[^;]+);?', val)
        )

    needs_close = not hasattr(fname, "close")
    fh: typing.TextIO = open(fname) if needs_close else fname
    try:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            seqname, source, feature, start, end, score, strand, frame, attr = (
                line.split("\t")
            )
            yield seqname, source, feature, int(start) - 1, int(end) - 1, float(
                "nan" if score == "." else score
            ), (None if strand == "." else strand), (
                -1 if frame == "." else int(frame)
            ), parse_attrs(
                attr
            )
    finally:
        if needs_close:
            fh.close()
