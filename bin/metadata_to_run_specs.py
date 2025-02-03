#!/usr/bin/env python

"""
Copyright © 2025 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
This file is part of DRAGoN.

This source code is licensed under the MIT License found in the
LICENSE file in the root directory of this source tree.
"""

import argparse
import os.path
import textwrap
import warnings

import numpy as np
import pandas as pd
from openpyxl.utils.cell import coordinate_to_tuple

# sys.tracebacklimit = 0


class DataInputError(Exception):
    pass


class CLI(argparse.Namespace):
    metadata: str
    sample_sheet: str = None

    _parser = argparse.ArgumentParser()
    _parser.add_argument(
        "metadata",
        help="Path to XLSX or JSON file describing the experiment collection",
    )
    _parser.add_argument(
        "-s",
        "--sample-sheet",
        dest="sample_sheet",
        help="Path to BCL2FASTQ sample sheet. "
        "If omitted, Excel sheet names or JSON keys must be prefix matches to the fastq files.",
    )

    _samplesheet_required_columns = [
        "Sample_ID",
        "Sample_Name",
        "Description",
        "index",
        "I7_Index_ID",
        "index2",
        "I5_Index_ID",
        "Sample_Project",
    ]
    _xlsx_required_columns = [
        "Plate",
        "PoolName",
        "Pool_Index_i7",
        "Well_Index",
        "WellBarcode",
        "SampleName",
        "Genome",
        "SampleID",
    ]
    _genome_species_map = {
        "Human": "GRCh38",
        "Mouse": "mm10",
        "Human + mouse": "GRCh38_mm10",
        "Macaque": "Mmul10",
        "Rat": "rat6",
    }
    _xlsx_required_unique_columns = ["Well_Index", "WellBarcode", "SampleID"]
    _xlsx_required_uniform_columns = ["Pool_Index_i7", "Genome"]
    _xlsx_future_required_columns = ["DataFolder"]

    def __init__(self, args=None):
        self.__class__._parser.parse_args(args, self)

    def parse_sample_sheet(self):
        if self.sample_sheet is None:
            return None
        sample_sheet = pd.read_csv(self.sample_sheet, header=None)
        if not (
            section_labels := sample_sheet.index[
                sample_sheet[0].str.startswith("[").replace(np.nan, False)
            ]
        ).any():
            raise DataInputError("Samplesheet missing section labels")
        section_starts = section_labels + 1
        section_ends = (section_labels[1:] - 1).append(pd.Index([len(sample_sheet)]))
        if not (data_mask := sample_sheet[0][section_labels] == "[Data]").any():
            raise DataInputError("Samplesheet missing required [Data] section")
        data_section = sample_sheet[
            section_starts[data_mask][0] + 1 : section_ends[data_mask][0]
        ].reset_index(drop=True)
        data_section.columns = sample_sheet.loc[section_starts[data_mask][0]]
        if missing_columns := [
            name
            for name in self._samplesheet_required_columns
            if name not in data_section
        ]:
            raise DataInputError(
                "Samplesheet missing required column(s) in [Data] section: "
                + ", ".join(missing_columns)
            )
        return data_section

    def parse_metadata_xlsx(self, sample_sheet=None):
        warned_species_deprecation = False
        input_error_group = []
        for sheet_name, sheet_data in pd.read_excel(
            self.metadata, None, header=2
        ).items():
            sheet_data.rename(
                columns=lambda name: name.strip().replace(" ", "_"), inplace=True
            )
            found_columns = set(self._xlsx_required_columns)
            if missing_columns := {
                name for name in found_columns if name not in sheet_data
            }:
                if "Genome" in missing_columns:
                    if not warned_species_deprecation:
                        warnings.warn(
                            'Column "Species" is deprecated and will be prohibited in a future version. '
                            'Use column "Genome" instead. '
                            "This warning will print once per experiment",
                            DeprecationWarning,
                        )
                        warned_species_deprecation = True
                    if unknown_species := set(sheet_data.Species) - set(
                        self._genome_species_map
                    ):
                        warnings.warn(
                            "Unrecognized species: " + ", ".join(unknown_species)
                        )
                    sheet_data["Genome"] = sheet_data.Species.apply(
                        lambda x: self._genome_species_map.get(x, f"Unknown__{x}")
                    )
                    missing_columns.remove("Genome")
                found_columns -= missing_columns
                if missing_columns:
                    input_error_group.append(
                        f"Metadata xlsx missing required column(s) in sheet {sheet_name}: "
                        + ", ".join(missing_columns)
                    )
                    continue  # fatal error
            for column in self._xlsx_future_required_columns:
                if column not in sheet_data:
                    warnings.warn(
                        'Column "DataFolder" will be required in a future version. This should be '
                        "the same as --IO.outdir."
                    )
            if empty_req_columns := [
                name for name in found_columns if sheet_data[name].isna().all()
            ]:
                input_error_group.append(
                    f"Metadata xlsx missing data in required column(s) in sheet {sheet_name}: "
                    + ", ".join(empty_req_columns)
                )
                continue  # fatal error
            sheet_data.dropna(inplace=True, subset=found_columns)
            if sheet_data.empty:
                input_error_group.append(
                    f"Metadata xlsx has no rows without missing values in required columns in sheet {sheet_name}"
                )
                continue  # fatal error
            if (
                illegal_chars := sheet_data[list(found_columns)]
                .select_dtypes(include="object")
                .drop("Pool_Index_i7", axis=1)
                .astype(str)
                .apply(lambda col: ~col.str.fullmatch(r"[\w.%-]+"))
            ).any(axis=None):
                wheres = textwrap.indent(
                    "\n".join(
                        illegal_chars[illegal_chars]
                        .stack()
                        .index.map(
                            lambda idx: 'Row {}, column {}: "{}"'.format(
                                *idx, sheet_data.loc[idx]
                            )
                        )
                    ),
                    "    ",
                )
                input_error_group.append(
                    f"Metadata xlsx contains values with illegal characters in sheet {sheet_name}. These may be invisible in the spreadsheet but show up in the formula bar\n{wheres}"
                )

            for column in self._xlsx_required_unique_columns:
                if column in sheet_data and sheet_data[column].nunique() != len(
                    sheet_data
                ):
                    input_error_group.append(
                        f"Duplicated values illegally found in sheet {sheet_name}, column {column}"
                    )

            for column in self._xlsx_required_uniform_columns:
                if column in sheet_data and sheet_data[column].nunique() != 1:
                    input_error_group.append(
                        f"Multiple values illegally found in sheet {sheet_name}, column {column}"
                    )

            sample_name = sheet_name
            i7_index = sheet_data["Pool_Index_i7"].iloc[0]
            if sample_sheet is not None:
                try:
                    sample_sheet_row = sample_sheet[
                        sample_sheet.I7_Index_ID == i7_index
                    ].iloc[0]
                except IndexError:
                    input_error_group.append(
                        f'Metadata sheet {sheet_name} specifies I7 index "{i7_index}" not found in Samplesheet'
                    )
                else:
                    sample_name = sample_sheet_row.Sample_Name

            if "MaskWells" in sheet_data:
                try:
                    if (
                        set(sheet_data["MaskWells"].dropna().astype(int).unique())
                        | {0, 1}
                    ) != {0, 1}:
                        raise TypeError(
                            f"incompatible dtype {sheet_data['MaskWells'].dtype}"
                        )
                    mask = sheet_data["MaskWells"].fillna(False).astype(bool)
                except TypeError as e:
                    input_error_group.append(
                        f"Metadata sheet {sheet_name} has MaskWells column but it could not be parsed as a true/false array: {e}"
                    )
            else:
                mask = pd.Series(False, index=sheet_data.index)

            try:
                barcodes = pd.concat(
                    [
                        sheet_data[["WellBarcode", "SampleID"]],
                        sheet_data["Well_Index"].apply(
                            lambda x: pd.Series(coordinate_to_tuple(x))
                        ),
                        mask,
                    ],
                    axis=1,
                ).dropna()
            except (ValueError, KeyError):
                input_error_group.append(
                    f"Failed to convert one or more Well_Index in sheet {sheet_name} to coordinates"
                )
                barcodes = None

            if input_error_group:
                continue

            barcodes_filename = f"{sample_name}_barcodes.txt"
            barcodes.to_csv(barcodes_filename, sep="\t", header=False, index=False)

        if input_error_group:
            raise DataInputError(
                "The following problems were detected in your input metadata:\n"
                + textwrap.indent("\n".join(input_error_group), "    ")
            )

    def parse_metadata_json(self, sample_sheet=None):
        # TODO: Define a spec and implement this parser
        raise NotImplementedError

    def main(self):
        sample_sheet = self.parse_sample_sheet()
        _, ext = os.path.splitext(self.metadata)
        if ext == ".xlsx":
            return self.parse_metadata_xlsx(sample_sheet=sample_sheet)
        elif ext == ".json":
            return self.parse_metadata_json(sample_sheet=sample_sheet)
        else:
            raise NotImplementedError(f"Unrecognized metadata extension: {ext}")


if __name__ == "__main__":
    CLI().main()
