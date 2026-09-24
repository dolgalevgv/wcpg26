#!/usr/bin/env python
"""Collect batch tables; apply gene-level FDR and retain all tested pairs."""

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control


GENE_FIELDS = [
    "n_variants", "ACAT_p", "p_adj_global", "p_adj_phenotype", "p_adj_gene",
    "is_egene", "fdr_scope", "fdr_threshold",
]


def read_genes(input_dir, threshold, scope):
    frames = []
    for path in sorted(input_dir.glob("*/gene_results.tsv")):
        frame = pd.read_csv(path, sep="\t", dtype={
            "phenotype_name": str, "gene_id": str, "top_MarkerID": str,
        })
        if frame.empty or frame["phenotype_name"].nunique() != 1:
            raise ValueError(f"Expected genes from one phenotype: {path}")
        frame["_batch_dir"] = str(path.parent)
        frames.append(frame)
    if not frames:
        raise ValueError(f"No batch gene summaries in {input_dir}")
    genes = pd.concat(frames, ignore_index=True)
    if genes.duplicated(["phenotype_name", "gene_id"]).any():
        raise ValueError("Duplicate phenotype/gene combinations across batches")
    if not genes["status"].isin(["tested", "no_testable_variants"]).all():
        raise ValueError("Unexpected gene status")
    genes["ACAT_p"] = pd.to_numeric(genes["ACAT_p"], errors="raise")
    genes["n_variants"] = pd.to_numeric(genes["n_variants"], errors="raise")
    tested = genes["status"].eq("tested")
    p = genes.loc[tested, "ACAT_p"]
    if not (np.isfinite(p) & p.between(0, 1)).all():
        raise ValueError("Tested genes require finite ACAT p-values in [0, 1]")
    if genes.loc[~tested, "ACAT_p"].notna().any():
        raise ValueError("Skipped genes unexpectedly have ACAT p-values")

    genes["p_adj_global"] = np.nan
    genes["p_adj_phenotype"] = np.nan
    if tested.any():
        genes.loc[tested, "p_adj_global"] = false_discovery_control(
            p.to_numpy(), method="bh"
        )
        for _, group in genes.loc[tested].groupby("phenotype_name"):
            genes.loc[group.index, "p_adj_phenotype"] = false_discovery_control(
                group["ACAT_p"].to_numpy(), method="bh"
            )
    genes["p_adj_gene"] = genes[f"p_adj_{scope}"]
    genes["is_egene"] = tested & genes["p_adj_gene"].le(threshold)
    genes["fdr_scope"] = scope
    genes["fdr_threshold"] = threshold
    return genes


def collect_pairs(genes, output_dir):
    association_dir = output_dir / "cis_associations"
    association_dir.mkdir(parents=True, exist_ok=True)
    required = {
        "phenotype_name", "gene_id", "CHR", "POS", "MarkerID",
        "Allele1", "Allele2", "BETA", "SE", "p.value",
    }
    pair_columns = None
    lead_frames = []
    for phenotype, group in genes.groupby("phenotype_name", sort=True):
        if "/" in phenotype or "\\" in phenotype:
            raise ValueError(f"Invalid phenotype filename: {phenotype}")
        with gzip.open(association_dir / f"{phenotype}.cis.tsv.gz", "wt") as out:
            write_header = True
            for batch_dir, batch_genes in group.groupby("_batch_dir", sort=True):
                path = Path(batch_dir) / "associations.tsv.gz"
                pairs = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
                if not required.issubset(pairs.columns):
                    raise ValueError(f"Missing association columns: {path}")
                if not pairs["phenotype_name"].eq(phenotype).all():
                    raise ValueError(f"Association phenotype mismatch: {path}")
                if not pairs["gene_id"].isin(batch_genes["gene_id"]).all():
                    raise ValueError(f"Associations without gene summaries: {path}")
                counts = pairs.groupby("gene_id").size().reindex(
                    batch_genes["gene_id"], fill_value=0
                )
                if not np.array_equal(counts.to_numpy(), batch_genes["n_variants"]):
                    raise ValueError(f"Association counts differ from summaries: {path}")
                if not np.array_equal(
                    counts.to_numpy() > 0, batch_genes["status"].eq("tested")
                ):
                    raise ValueError(f"Association counts disagree with status: {path}")

                p = pd.to_numeric(
                    pairs["p.value"].replace({"NA": np.nan, "": np.nan}),
                    errors="raise",
                )
                if (p.notna() & ~p.between(0, 1)).any():
                    raise ValueError(f"Invalid association p-values: {path}")
                tested_genes = set(batch_genes.loc[
                    batch_genes["status"].eq("tested"), "gene_id"
                ])
                if tested_genes - set(pairs.loc[p.notna(), "gene_id"]):
                    raise ValueError(f"Tested genes without valid p-values: {path}")

                pairs.insert(2, "variant_id", (
                    pairs["CHR"] + ":" + pairs["POS"] + ":"
                    + pairs["Allele1"] + ":" + pairs["Allele2"]
                ))
                if pair_columns is None:
                    pair_columns = list(pairs.columns)
                elif list(pairs.columns) != pair_columns:
                    raise ValueError(f"Association columns differ between batches: {path}")
                pairs.to_csv(out, sep="\t", index=False, header=write_header)
                write_header = False

                significant = batch_genes.loc[batch_genes["is_egene"]]
                select = p.notna() & pairs["gene_id"].isin(significant["gene_id"])
                if select.any():
                    indices = p[select].groupby(pairs.loc[select, "gene_id"]).idxmin()
                    lead_frames.append(pairs.loc[indices].merge(
                        significant[["phenotype_name", "gene_id", *GENE_FIELDS]],
                        on=["phenotype_name", "gene_id"], validate="one_to_one",
                    ))
    leads = (pd.concat(lead_frames, ignore_index=True) if lead_frames else
             pd.DataFrame(columns=[*pair_columns, *GENE_FIELDS]))
    if len(leads) != int(genes["is_egene"].sum()):
        raise ValueError("Missing lead associations for significant genes")
    leads.to_csv(output_dir / "lead_eqtl.tsv", sep="\t", index=False, na_rep="NA")


def collect_results(input_dir, threshold=0.05, scope="global", output_dir=Path(".")):
    if not 0 < threshold < 1 or scope not in {"global", "phenotype"}:
        raise ValueError("Invalid FDR threshold or scope")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    genes = read_genes(Path(input_dir), threshold, scope)
    collect_pairs(genes, output_dir)
    genes = genes.drop(columns="_batch_dir").sort_values(["phenotype_name", "gene_id"])
    genes.to_csv(output_dir / "all_genes.tsv", sep="\t", index=False, na_rep="NA")
    genes.loc[genes["is_egene"]].to_csv(
        output_dir / "significant_genes.tsv", sep="\t", index=False, na_rep="NA"
    )
    summary = genes.groupby("phenotype_name", sort=True).agg(
        n_genes=("gene_id", "size"), n_tested=("ACAT_p", "count"),
        n_pairs=("n_variants", "sum"), n_egenes=("is_egene", "sum"),
    )
    summary["n_no_testable_variants"] = summary["n_genes"] - summary["n_tested"]
    summary.reset_index().to_csv(output_dir / "phenotype_summary.tsv", sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("fdr_threshold", type=float)
    parser.add_argument("--fdr-scope", choices=["global", "phenotype"], default="global")
    args = parser.parse_args()
    collect_results(args.input_dir, args.fdr_threshold, args.fdr_scope)


if __name__ == "__main__":
    main()

