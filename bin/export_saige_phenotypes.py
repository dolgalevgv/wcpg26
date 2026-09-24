#!/usr/bin/env python
"""Export deterministic gene batches, retaining one input archive per batch."""

import argparse
import re
import tarfile
import tempfile
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import sparse


def export_batches(adata_path, donors_path, batch_size, output_dir=Path(".")):
    if batch_size < 1:
        raise ValueError("Batch size must be positive")

    output_dir = Path(output_dir)
    batch_dir = output_dir / "batches"
    batch_dir.mkdir(parents=True, exist_ok=True)
    donors = pd.read_csv(donors_path).drop(columns=["donor", "n_cells"])
    adata = ad.read_h5ad(adata_path)
    if not adata.var_names.is_unique or adata.n_vars == 0:
        raise ValueError("Expected at least one gene and unique gene IDs")
    if any(not re.fullmatch(r"[A-Za-z0-9_.-]+", gene)
           or gene in {".", ".."} for gene in adata.var_names):
        raise ValueError("Gene IDs must be safe filenames")

    order = sorted(range(adata.n_vars), key=lambda i: adata.var_names[i])
    regions = adata.var[["chrom", "cis_start", "cis_end"]].copy()
    regions.index.name = "gene_id"
    regions["phenotype_file"] = regions.index.astype(str) + ".tsv.gz"
    if sparse.issparse(adata.X):
        adata.X = adata.X.tocsc()
    expression_matrix = adata.X
    batches = []

    # Only the current gene TSV exists outside the archives. This also keeps
    # the persistent file count low when node-local scratch is unavailable.
    with tempfile.TemporaryDirectory(prefix="export-", dir=output_dir) as tmp:
        tmp = Path(tmp)
        covariates = tmp / "covariates.txt"
        covariates.write_text(",".join(donors.columns) + "\n")
        for start in range(0, len(order), batch_size):
            indices = order[start:start + batch_size]
            batch_id = f"batch{start // batch_size + 1:04d}"
            archive_name = f"{batch_id}.tar"
            manifest = tmp / "manifest.tsv"
            regions.iloc[indices].reset_index().to_csv(
                manifest, sep="\t", index=False
            )
            with tarfile.open(batch_dir / archive_name, "w") as archive:
                archive.add(manifest, arcname="manifest.tsv")
                archive.add(covariates, arcname="covariates.txt")
                for i in indices:
                    gene = adata.var_names[i]
                    column = expression_matrix[:, i]
                    expression = (column.toarray().ravel()
                                  if sparse.issparse(column) else column.ravel())
                    table = adata.obs.assign(expression=expression).reset_index(
                        names="cell_id"
                    )
                    phenotype_file = tmp / f"{gene}.tsv.gz"
                    table.to_csv(phenotype_file, sep="\t", index=False)
                    archive.add(phenotype_file, arcname=phenotype_file.name)
                    phenotype_file.unlink()
            batches.append({
                "batch_id": batch_id,
                "archive": archive_name,
                "n_genes": len(indices),
            })
    pd.DataFrame(batches).to_csv(output_dir / "manifest.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("adata", type=Path)
    parser.add_argument("donors", type=Path)
    parser.add_argument("--batch-size", type=int, default=100)
    args = parser.parse_args()
    export_batches(args.adata, args.donors, args.batch_size)


if __name__ == "__main__":
    main()
