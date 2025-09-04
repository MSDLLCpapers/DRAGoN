#!/usr/bin/env python

"""
Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
This file is part of DRAGoN.

This source code is licensed under the MIT License found in the
LICENSE file in the root directory of this source tree.
"""

import abc
import argparse
import collections.abc
import contextlib
import functools
import json
import logging
import multiprocessing
import multiprocessing.sharedctypes
import multiprocessing.synchronize
import pathlib
import sys
import textwrap
import time
import typing
from typing import Literal, TypedDict

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn

try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None


UNMATCHED_REFRESH_FREQUENCY = 1 / 60

logging.basicConfig(
    level=logging.INFO, format="[%(asctime)s] %(name)s: %(levelname)s: %(message)s"
)
log = logging.getLogger(pathlib.Path(__file__).stem)


@contextlib.contextmanager
def figure(*args, **kwargs):
    fig = plt.figure(*args, **kwargs)
    yield fig
    plt.close(fig)


class DemuxReport(TypedDict):
    read_count: int
    well_counts: dict[str, int]
    written_counts: list[int]
    mismatch_counts: list[dict[str, int]]
    fail_counts: dict[str, dict[str, int]]
    unmatched_counts: dict[str, int]
    plate_map: list[list[int]]
    mapped_wells: list[list[bool]]


def load_demux_report(filenames: list[str]) -> DemuxReport:
    """
    Loads and consolidates the demultiplex reports
    :param filenames: List of filenames to the demultiplex json statistics
    :return: A 2-tuple containing a Series of total counts per barcode and a dict summing up all the input json reports. The dict is also saved to 'demux.json'.
    """
    with contextlib.ExitStack() as stack:
        json_orig: list[DemuxReport] = [
            json.load(stack.enter_context(open(fname))) for fname in filenames
        ]
    demux_consolidated: DemuxReport = {
        "read_count": sum(x["read_count"] for x in json_orig),
        "well_counts": pd.DataFrame([x["well_counts"] for x in json_orig])
        .sum()
        .to_dict(),
        "written_counts": pd.DataFrame([x["written_counts"] for x in json_orig])
        .sum()
        .to_list(),
        "mismatch_counts": [
            pd.DataFrame([x["mismatch_counts"][i] for x in json_orig]).sum().to_dict()
            for i in range(len(json_orig[0]["mismatch_counts"]))
        ],
        "fail_counts": {
            i: pd.DataFrame([x["fail_counts"][i] for x in json_orig]).sum().to_dict()
            for i in json_orig[0]["fail_counts"]
        },
        "unmatched_counts": functools.reduce(
            lambda a, b: a.add(b, fill_value=0),
            (pd.Series(x["unmatched_counts"]) for x in json_orig),
            pd.Series(dtype=int),
        )
        .sort_values(ascending=False)
        .astype(int)
        .to_dict(),
    }
    return demux_consolidated


def read_barcodes(filename: str) -> pd.Series:
    df = pd.read_csv(
        filename,
        sep="\t",
        names=["seq", "name", "x", "y", "mask"],
        index_col="name",
    )
    df["mask"] = df["mask"].astype(bool)
    return df


class HeuristicMeta(type):
    _children: "set[type[Heuristic]]" = set()

    def __iter__(cls):
        return iter(cls._children)


class Heuristic(metaclass=HeuristicMeta):
    def __init_subclass__(cls, **kwargs):
        Heuristic._children.add(cls)

    def __init__(self, demux_report: DemuxReport, mask: pd.Series):
        self.demux_report = demux_report
        self.well_counts = pd.Series(demux_report["well_counts"])
        self.fail_counts = pd.DataFrame(demux_report["fail_counts"])
        self.failed: pd.Series = None
        self.mask = mask.copy()
        self.mask.loc["noBC"] = False

    @property
    @abc.abstractmethod
    def limit(self) -> float: ...

    @property
    @abc.abstractmethod
    def message(self) -> str: ...

    def __call__(self):
        self.failed = self.well_counts.loc[~self.mask] > self.limit
        return self.failed.any()

    def inv(self, count: int) -> str:
        return "check demux.json"

    def __str__(self):
        if not self.failed.any():
            return ""
        join = textwrap.indent(
            "\n".join(
                f"{name} = {count:,d} ({self.inv(count):s})"
                for name, count in self.well_counts.loc[~self.mask]
                .loc[self.failed]
                .items()
            ),
            " " * 4,
        )
        return f"{self.message}:\n{join}"


class IQR(Heuristic):
    def __init__(self, demux_report: DemuxReport, mask: pd.Series):
        super().__init__(demux_report, mask)
        self.well_counts.drop("noBC", inplace=True)
        self.mask.drop("noBC", inplace=True)

    @functools.cached_property
    def lower_quantile(self):
        return self.well_counts.loc[~self.mask].quantile(0.25)

    @functools.cached_property
    def upper_quantile(self):
        return self.well_counts.loc[~self.mask].quantile(0.75)

    @property
    def limit(self):
        return self.upper_quantile + 1.75 * (self.upper_quantile - self.lower_quantile)

    @property
    def message(self):
        return f"The following wells have excessive read counts (> 75th pctl + 1.75 * iqr = {self.limit:,.2f})"

    def inv(self, count: int):
        percent_over = (
            100
            * (count - self.upper_quantile)
            / (self.upper_quantile - self.lower_quantile)
        )
        return f"{percent_over:,.2f}% of IQR over upper quantile = {self.upper_quantile:,.2f}"


class Fraction(Heuristic):
    @property
    def limit(self):
        return 0.15 * self.well_counts.loc[~self.mask].sum()

    @property
    def message(self):
        return f"The following wells have excessive read counts (> 15% experiment size = {self.limit:,.2f})"

    def inv(self, count: int):
        pct = 100 * count / self.well_counts.loc[~self.mask].sum()
        return f"{pct:,.2f}% of total reads"


class TooManyFailed(Heuristic):
    @property
    def limit(self):
        return 0.5 * self.well_counts.loc[~self.mask]

    def __call__(self):
        self.failed = self.fail_counts.loc[~self.mask, "PASS"] < self.limit
        return self.failed.any()

    @property
    def message(self):
        return "The following wells had more than 50%% of reads fail QC"


def mp_init(
    mapped2coords: dict[str, list[int, int]],
    mismatch: int,
    queue: multiprocessing.Queue,
):
    mp_init.mapped2coords = list(mapped2coords)
    mp_init.mismatch = mismatch
    mp_init.queue = queue


def match_barcodes(seqA: str, mapped2coords: collections.abc.Container, mismatch: int):
    cur = {candidate: 0 for candidate in mapped2coords if candidate}
    for i, c in enumerate(seqA):
        for candidate in list(cur):
            if candidate[i] != c:
                cur[candidate] += 1
                if cur[candidate] > mismatch:
                    cur.pop(candidate)
        if not cur:
            return None, None
    return min(cur.items(), key=lambda t: t[1])


def match_barcodes_worker(barcode: str, count: int):
    candidate, distance = match_barcodes(
        barcode, mp_init.mapped2coords, mp_init.mismatch
    )
    mp_init.queue.put((barcode, count, candidate, distance))


class CLI(argparse.Namespace):
    dmux: list[str]
    layout: typing.TextIO | None
    handle_outlier_wells: Literal["ignore", "warn", "raise"] = "raise"
    mismatch: int = 1
    xoff: int = 0
    yoff: int = 0
    barcodes: typing.TextIO
    dont_resolve_nobc: bool = False
    plot_masked_wells: bool = True
    nobc_resolve_threads: int = 1

    def __init__(self, args=None):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "dmux",
            nargs="+",
            help="JSON reports from demultiplex.py (Preprocessing:Demultiplex)",
        )
        parser.add_argument("--layout", type=argparse.FileType())
        parser.add_argument(
            "--handle-outlier-wells",
            choices=["ignore", "warn", "raise"],
            default=CLI.handle_outlier_wells,
        )
        parser.add_argument("--xoff", type=int, default=CLI.xoff)
        parser.add_argument("--yoff", type=int, default=CLI.yoff)
        parser.add_argument("--mismatch", type=int, default=CLI.mismatch)
        parser.add_argument("--barcodes", type=argparse.FileType(), required=True)
        parser.add_argument("--dont-resolve-nobc", action="store_true")
        parser.add_argument(
            "--dont-plot-masked-wells", action="store_false", dest="plot_masked_wells"
        )
        parser.add_argument(
            "--nobc-resolve-threads", type=int, default=CLI.nobc_resolve_threads
        )
        _, self.remaining_args = parser.parse_known_args(args, self)
        if self.remaining_args:
            log.warn(
                "extra commandline arguments not recognized: %s",
                " ".join(self.remaining_args),
            )

    def make_grid(self, demux_report: DemuxReport, barcode_coords: pd.DataFrame):
        if self.layout is None:
            return

        layout = [row.rstrip("\n").split(",") for row in self.layout]
        nrow = len(layout)
        ncol = layout and len(layout[0])

        mapped_counts: list[list[int]] = [[0 for _ in range(ncol)] for _ in range(nrow)]
        mapped_wells: list[list[str | None]] = [
            [None for _ in range(ncol)] for _ in range(nrow)
        ]
        metadata_barcodes: set[tuple[int, int]] = set()
        for key, value in demux_report["well_counts"].items():
            if key not in barcode_coords.index:
                continue
            x, y = (barcode_coords.loc[key] - 1).to_list()
            if 0 <= y < nrow and 0 <= x < ncol:
                mapped_counts[y][x] += value
                metadata_barcodes.add((y, x))
                mapped_wells[y][x] = key
        mapped2coords = {
            barcode: (y, x)
            for y, row in enumerate(layout)
            for x, barcode in enumerate(row)
            if (y, x) not in metadata_barcodes
        }
        for barcode, (y, x) in mapped2coords.items():
            mapped_counts[y][x] += demux_report["unmatched_counts"].pop(barcode, 0)

        n_unmatched = len(demux_report["unmatched_counts"])
        print(n_unmatched, "unmatched barcodes", file=sys.stderr)

        if not self.dont_resolve_nobc:
            is_atty = sys.stdout.isatty()
            start_time = time.monotonic()
            next_update = multiprocessing.Value("d", start_time)

            def on_starmap_return(
                i: int,
                barcode: str,
                count: int,
                candidate: str | None,
                distance: int | None,
            ):
                if is_atty:
                    with next_update:
                        cur_time = time.monotonic()
                        if i in {1, n_unmatched} or cur_time >= next_update.value:
                            print(f"\r{i}\t{barcode}", end="", file=sys.stderr)
                            while cur_time >= next_update.value:
                                next_update.value += UNMATCHED_REFRESH_FREQUENCY

                if distance is not None:
                    y, x = mapped2coords[candidate]
                    mapped_counts[y][x] += count
                    demux_report["unmatched_counts"].pop(barcode)

            if self.nobc_resolve_threads > 1:
                queue = multiprocessing.Queue()
                with multiprocessing.Pool(
                    processes=self.nobc_resolve_threads,
                    initializer=mp_init,
                    initargs=(mapped2coords, self.mismatch, queue),
                ) as pool:
                    result = pool.starmap_async(
                        match_barcodes_worker,
                        demux_report["unmatched_counts"].items(),
                    )
                    pool.close()
                    i = 1
                    while not result.ready() or not queue.empty():
                        barcode, count, candidate, distance = queue.get()
                        on_starmap_return(i, barcode, count, candidate, distance)
                        i += 1
                    queue.close()
                    pool.join()
                    queue.join_thread()
            else:
                for i, (barcode, count) in enumerate(
                    demux_report["unmatched_counts"].copy().items(),
                    1,
                ):
                    candidate, distance = match_barcodes(
                        barcode, mapped2coords, self.mismatch
                    )
                    on_starmap_return(i, barcode, count, candidate, distance)
            print("", file=sys.stderr)

        # CI: Force consistent order
        demux_report["unmatched_counts"] = dict(
            sorted(demux_report["unmatched_counts"].items(), key=lambda t: t[::-1])
        )
        if not all(x is None for y in mapped_wells for x in y):
            demux_report["plate_map"] = mapped_counts
            demux_report["mapped_wells"] = mapped_wells
        else:
            logging.warn("Plate layout does not match experiment. Skipping plotting")
            demux_report["plate_map"] = None
            demux_report["mapped_wells"] = None

    def handle_errors(self, demux_report, mask):
        message = []
        for cls in Heuristic:
            heuristic = cls(demux_report, mask)
            if heuristic() and self.handle_outlier_wells != "ignore":
                message.append(str(heuristic))
            failed = heuristic.failed.to_frame("flag").reset_index(names="well_name")
            failed = failed.loc[failed["flag"], ["well_name"]]
            if isinstance(heuristic.limit, float):
                failed = pd.concat(
                    [
                        pd.DataFrame(
                            [str(round(heuristic.limit))],
                            columns=["well_name"],
                            index=["limit"],
                        ),
                        failed,
                    ],
                    axis=0,
                )
            failed.to_csv(f"{cls.__name__}.failed.txt", header=False, sep="\t")
        if message:
            print(
                "One or more data quality issues were detected (see below)\n\n"
                + textwrap.indent("\n\n".join(message), " " * 4)
            )

    def write_demux_report(self, demux_report: DemuxReport):
        with open("demux.json", "w") as ofp:
            json.dump(demux_report, ofp, indent=4)

    def platemap_to_csv(self, demux_report: DemuxReport):
        if self.layout is None:
            return

        mapped_counts = demux_report["plate_map"]
        with open("demux.layout.csv", "w") as ofp:
            print("", *range(1, 25), sep=",", file=ofp)
            for y, row in enumerate(mapped_counts):
                print(chr(y + ord("A")), *row, sep=",", file=ofp)

    def graphical_report(self, demux_report: DemuxReport, mask: pd.Series):
        if self.layout is None:
            return

        matplotlib.rc("axes", titlesize=24, labelsize=18)
        noBC = sum(demux_report["unmatched_counts"].values())
        plate_map = pd.DataFrame(
            demux_report["plate_map"],
            columns=range(1, 25),
            index=list("ABCDEFGHIJKLMNOP"),
        )
        wells_to_name = pd.DataFrame(
            demux_report["mapped_wells"],
            columns=range(1, 25),
            index=list("ABCDEFGHIJKLMNOP"),
            dtype=str,
        )
        long = pd.merge_ordered(
            wells_to_name.melt(var_name="column", value_name="name").reset_index(
                names="row"
            ),
            plate_map.melt(var_name="column", value_name="count").reset_index(
                names="row"
            ),
            on=["row", "column"],
            how="outer",
        )
        long["mask"] = (
            long.reset_index(names="index")
            .set_index("name")
            .merge(mask, left_index=True, right_index=True, how="outer")
            .set_index("index")["mask"]
            .astype(bool)
        )
        long["name"] = long["name"].fillna(
            long.apply(lambda row: f'UNASSIGNED__{row["row"]}{row["column"]}', axis=1),
        )
        long = pd.concat(
            [
                long,
                long["name"].str.extract(
                    r"^(?P<pool>[^_]+)_(?P<coords>[^_]+)__(?P<label>.+)$"
                ),
            ],
            axis=1,
        )
        wells_to_name = wells_to_name.where(wells_to_name.isin(mask.loc[~mask].index))
        if not self.plot_masked_wells:
            long = long.loc[~long["mask"]]
        vmax = plate_map.max(axis=None)
        if noBC > vmax:
            vmax *= 1.1
        plate_map.loc["-", 0] = noBC
        wells_to_name.loc["-", 0] = "noBC"

        # Heatmap
        with figure(figsize=[20, 12]) as fig:
            ax = fig.gca()
            seaborn.heatmap(
                plate_map,
                ax=ax,
                vmin=0,
                vmax=vmax,
                mask=wells_to_name.isna(),
                cbar=True,
                cbar_kws={"format": "%d"},
                cmap="Reds",
                linewidths=0.5,
                linecolor="0.69",
            )
            seaborn.heatmap(
                plate_map,
                ax=ax,
                vmin=0,
                vmax=vmax,
                mask=~wells_to_name.isna(),
                cbar=False,
                cmap="Greys",
                linewidths=0.5,
                linecolor="0.69",
            )
            ax.set_title("Demultiplexed plate map")
            ax.tick_params(labelsize=18)
            fig.savefig("plate_heatmap.png")

        # Violin plot
        with figure(figsize=[12, 20]) as fig:
            ax = fig.gca()
            seaborn.violinplot(long, y="count", inner="quart", ax=ax)
            q25, q75 = long.loc[~long["mask"], "count"].quantile((0.25, 0.75))
            reasonable_max = q75 + 1.75 * (q75 - q25)
            if (overfulls := long["count"] > reasonable_max).any():
                texts = long.loc[overfulls].apply(
                    lambda row: ax.text(
                        0,
                        row["count"],
                        "{label:s}\n({pool:s}:{coords:s})".format(
                            **row.astype("str").to_dict()
                        ),
                        ha="center",
                    ),
                    axis=1,
                )
                if adjust_text is not None:
                    adjust_text(
                        texts.to_list(),
                        arrowprops={"arrowstyle": "-", "color": "black", "lw": 0.5},
                    )
            ax.set_title("Reads per well")
            ax.tick_params(labelsize=18)
            ax.ticklabel_format(style="plain", axis="y")
            fig.savefig("violin_plot.png")

    def main(self):
        demux_report = load_demux_report(self.dmux)
        barcodes = read_barcodes(self.barcodes)
        self.make_grid(demux_report, barcodes[["x", "y"]])
        self.handle_errors(demux_report, barcodes["mask"])
        self.write_demux_report(demux_report)
        if demux_report["plate_map"] is not None:
            self.platemap_to_csv(demux_report)
            self.graphical_report(demux_report, barcodes["mask"])


if __name__ == "__main__":
    CLI().main()
