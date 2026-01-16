#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title:          format_coverage.py
Description:    Aggregates coverage statistics from TSV files into a multi-sheet Excel report.
                Calculates observed vs. expected percentages per amplicon.
"""

# %%
# from sys import argv
import os
import pandas as pd

PATH = "/home/james/repos/github/pfal_variant_calling/output/coverage/QC/results/"

# %%
reads = pd.read_csv(
    os.path.join(PATH, "readcount.tsv"),
    sep="\t",
)

readcount = pd.pivot_table(
    reads,
    values="reads",
    columns=["sample_id"],
    index=["amplicon"],
)

readcount["total"] = readcount.sum(axis=1)
readcount["obs % per amplicon"] = round(
    readcount["total"] / readcount["total"].sum() * 100, 2
)
readcount["exp % / amplicon"] = round(1 / readcount.shape[0] * 100)
readcount["fold diff"] = round(
    readcount["obs % per amplicon"] / readcount["exp % / amplicon"], 2
)

# %%

coverage = pd.read_csv(
    os.path.join(PATH, "coverage.tsv"),
    sep="\t",
)

coverage = pd.pivot_table(
    coverage,
    values="mean",
    columns=["sample_id"],
    index=["amplicon"],
)
coverage = coverage[coverage.index != "total"]
coverage.fillna(0, inplace=True)

# %%
basecount = pd.read_csv(
    os.path.join(PATH, "basecount.tsv"),
    sep="\t",
)

basecount = pd.pivot_table(
    basecount,
    values="nbases",
    columns=["sample_id"],
    index=["amplicon"],
)

# %%

with pd.ExcelWriter(
    os.path.join(PATH, "coverage_report.xlsx"),
    engine="openpyxl",
    mode="w",
) as writer:
    coverage.to_excel(writer, sheet_name="Mean coverage")
    readcount.to_excel(writer, sheet_name="Number of reads")
    basecount.to_excel(writer, sheet_name="Number of bases")

# %%
