#!/usr/bin/env python

from pathlib import Path
import sys

import anndata as ad
import pandas as pd


adata_path, donors_path = [Path(p) for p in sys.argv[1:3]]

donors = pd.read_csv(donors_path).drop(columns=["donor", "n_cells"])
covariates = donors.columns

with open("covariates.txt", "w") as f:
    f.write(",".join(covariates) + "\n")

adata = ad.io.read_h5ad(adata_path)

manifest = adata.var[["chrom", "cis_start", "cis_end"]].copy().reset_index(names="gene_id")
manifest["phenotype_file"] = manifest["gene_id"].astype(str) + ".tsv.gz"

manifest.to_csv("manifest.csv", index=False)

adata.X = adata.X.tocsc()

Path("gene_inputs").mkdir(exist_ok=True)

for i, gene in enumerate(adata.var.index):
    expression = adata.X[:, i].toarray().ravel()
    sub_gene = adata.obs.assign(expression=expression).reset_index(names="cell_id")
    sub_gene.to_csv(f"gene_inputs/{gene}.tsv.gz", sep="\t", index=False)
