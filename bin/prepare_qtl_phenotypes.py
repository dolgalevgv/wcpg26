#!/usr/bin/env python

import argparse
import itertools
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("adata_path")
parser.add_argument("donors_path")
parser.add_argument("pca_path")
parser.add_argument("regions_path")

parser.add_argument("--cell_type_col", type=str, required=True)
parser.add_argument("--treat_col", type=str, required=True)
parser.add_argument("--qcovar", type=str, required=True)
parser.add_argument("--min-n-cells", type=int, required=True)
parser.add_argument("--min-n-donors", type=int, required=True)
parser.add_argument("--min-cells-frac", type=float, required=True)
parser.add_argument("--min-donor-frac", type=float, required=True)

args = parser.parse_args()

adata_path = args.adata_path
donors_path = args.donors_path
pca_path = args.pca_path
regions_path = args.regions_path

cell_type_col = args.cell_type_col
treat_col = args.treat_col
qcovar = args.qcovar.split(",") if "," in args.qcovar else [args.qcovar]

min_n_donors = args.min_n_donors
min_n_cells = args.min_n_cells
min_cells_frac = args.min_cells_frac
min_donor_frac = args.min_donor_frac

adata = ad.io.read_h5ad(adata_path)
adata.obs["log_total_counts"] = np.log(adata.obs["total_counts"])
cell_type_levels = adata.obs[cell_type_col].unique()
treat_levels = adata.obs[treat_col].unique()

donors = pd.read_csv(donors_path, index_col=0)
pca = pd.read_csv(pca_path, sep="\t", index_col=0)
regions = pd.read_csv(regions_path, sep="\t", index_col=0)

donors = donors.drop(columns="vcf_id").join(pca)

for c in qcovar:
    onehot = pd.get_dummies(
        donors[c],
        prefix=c,
        drop_first=True,
        dtype=int,
    )

    donors = donors.drop(columns=c).join(onehot)

manifest = []

for t, c in itertools.product(treat_levels, cell_type_levels):
    sub_adata = adata[adata.obs[treat_col].eq(t) & adata.obs[cell_type_col].eq(c)]
    sub_adata.obs = sub_adata.obs[["donor", "cell_type", "treat", "log_total_counts"]]

    n_cells = sub_adata.obs.groupby("donor").size().to_frame("n_cells")
    sub_donors = donors.join(n_cells)
    sub_donors["n_cells"] = sub_donors["n_cells"].fillna(0).astype(int)

    donors_excluded = sub_donors.loc[sub_donors["n_cells"] < min_n_cells]
    sub_donors = sub_donors.loc[sub_donors["n_cells"] >= min_n_cells]

    if sub_donors.shape[0] < min_n_donors:
        print(f"Too few donors for {c}_{t} comparison: {sub_donors.shape[0]} < {min_n_donors}, skipping it")
        continue

    sub_adata = sub_adata[sub_adata.obs["donor"].isin(sub_donors.index)]

    sub_adata.obs = sub_adata.obs.join(
        sub_donors, 
        on="donor",
        how="left",
        validate="many_to_one",
    )
    sub_adata.obs = sub_adata.obs.drop(columns=["n_cells", cell_type_col, treat_col])

    sub_adata.var = sub_adata.var.set_index("gene_id").loc[:, ["mt"]]
    regions_keep = sub_adata.var.index.isin(regions.index)
    mt_keep = ~sub_adata.var["mt"]
    
    sub_counts = sub_adata.layers["counts"].tocsr()

    donor_frac = {}

    for donor in sub_donors.index:
        cells_donor = (sub_adata.obs["donor"] == donor).to_numpy()
        sub_counts_donor = sub_counts[cells_donor]
        donor_frac[donor] = sub_counts_donor.count_nonzero(axis=0) / sub_counts_donor.shape[0]

    donor_frac = pd.DataFrame(donor_frac, index=sub_adata.var.index)
    donor_frac_genes = (donor_frac >= min_cells_frac).sum(axis=1) / len(sub_donors.index)
    donor_keep = (donor_frac_genes >= min_donor_frac)

    genes = sub_adata.var.copy()
    sub_adata.var = sub_adata.var.drop(columns="mt")
    genes["in_regions"] = regions_keep
    genes["donor_frac"] = donor_frac_genes
    genes["donor_frac_keep"] = donor_keep

    gene_keep = regions_keep & mt_keep & donor_keep
    genes_excluded = genes.loc[~gene_keep]

    sub_adata = sub_adata[:, gene_keep]
    sub_adata.var = sub_adata.var.join(regions)

    manifest.append([c, t, sub_donors.shape[0], sub_adata.shape[1]])

    stratum = Path(f"{c}_{t}")
    Path(stratum).mkdir(exist_ok=True)

    # Rebuild to simplify the object
    sub_adata = ad.AnnData(
        X=sub_adata.layers["counts"].copy(),
        obs=sub_adata.obs.copy(),
        var=sub_adata.var.copy()
    )

    ad.io.write_h5ad(stratum / Path("phenotypes.h5ad"), sub_adata)

    sub_donors.reset_index(names="donor").to_csv(stratum / Path("donors.csv"), index=False)
    if not donors_excluded.empty:
        donors_excluded.reset_index(names="donor").to_csv(stratum / Path("donors_excluded.csv"), index=False)

    genes_excluded.reset_index(names="gene_id").to_csv(stratum / Path("genes_excluded.csv"), index=False)
    sub_adata.var.reset_index(names="gene_id").to_csv(stratum / Path("genes.csv"), index=False)

manifest = pd.DataFrame(manifest, columns=[cell_type_col, treat_col, "n_donors", "n_genes"])
manifest.to_csv("manifest.csv", index=False)