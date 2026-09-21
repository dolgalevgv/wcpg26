#!/usr/bin/env python

import sys
from pathlib import Path

import numpy as np
import pandas as pd


gtf_path, fai_path, regions_path = [Path(path) for path in sys.argv[1:4]]
cis_window = int(sys.argv[4]) if len(sys.argv) > 4 else 1000000

chr_lengths = pd.read_csv(fai_path, sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"], dtype={"chrom": str})
chr_lengths = chr_lengths.loc[chr_lengths["chrom"].isin([str(i) for i in range(1, 23)])]
chr_lengths = dict(zip(chr_lengths["chrom"].to_list(), chr_lengths["length"].to_list()))

gtf = pd.read_csv(gtf_path, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8], names=["chrom", "type", "start", "end", "strand", "attrs"], dtype={"chrom": str})
gtf = gtf.loc[gtf["type"] == "gene"]
gtf = gtf.loc[gtf["chrom"].isin(chr_lengths)]
gtf = gtf.loc[gtf["strand"].isin(["+", "-"])]

attrs = (
    gtf["attrs"]
    .str.findall(r'(\w+)\s+"([^"]*)"')
    .map(dict)
    .apply(pd.Series)
)
gtf = gtf.drop(columns="attrs").join(attrs).loc[:, ["gene_id", "gene_name", "chrom", "strand", "start", "end"]]

gtf["tss"] = np.where(gtf["strand"] == "+", gtf["start"], gtf["end"])

groups = []
for chrom, group in gtf.groupby("chrom"):
    group["cis_start"] = (group["tss"] - cis_window).clip(lower=1).astype(int)
    group["cis_end"] = (group["tss"] + cis_window).clip(upper=chr_lengths[chrom]).astype(int)

    groups.append(group)

gene_regions = pd.concat(groups, ignore_index=True)

gene_regions.to_csv(regions_path, sep="\t", index=False, header=True)