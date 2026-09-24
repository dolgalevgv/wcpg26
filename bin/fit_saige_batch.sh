#!/usr/bin/env bash
# Arguments are absolute paths, followed by the existing SAIGE fit settings.
set -euo pipefail
inputs=$1
bed=$2
fam=$3
bim=$4
output=$5
tol=$6
maxiter=$7
threads=$8

batch_tmp=$(mktemp -d "${PWD}/saige-fit.XXXXXX")
trap 'rm -rf -- "$batch_tmp"' EXIT
tar -xf "$inputs" -C "$batch_tmp"
mkdir "$batch_tmp/models"
covariate_columns=$(cat "$batch_tmp/covariates.txt")
tar -cf "$output" -C "$batch_tmp" manifest.tsv

tail -n +2 "$batch_tmp/manifest.tsv" |
while IFS=$'\t' read -r gene chrom cis_start cis_end phenotype_file; do
    printf 'Fitting %s\n' "$gene"
    prefix="$batch_tmp/models/$gene"
    step1_fitNULLGLMM_qtl.R \
        --phenoFile="$batch_tmp/$phenotype_file" \
        --phenoCol=expression \
        --sampleIDColinphenoFile=donor \
        --cellIDColinphenoFile=cell_id \
        --covarColList="$covariate_columns" \
        --sampleCovarColList="$covariate_columns" \
        --offsetCol=log_total_counts \
        --traitType=count \
        --bedFile="$bed" --famFile="$fam" --bimFile="$bim" \
        --useGRMtoFitNULL=FALSE --useSparseGRMtoFitNULL=FALSE \
        --LOCO=FALSE --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE --isCovariateTransform=TRUE \
        --skipModelFitting=FALSE --skipVarianceRatioEstimation=FALSE \
        --isCateVarianceRatio=FALSE --IsOverwriteVarianceRatioFile=TRUE \
        --isStoreSigma=TRUE --isShrinkModelOutput=TRUE \
        --tol="$tol" --maxiter="$maxiter" --nThreads="$threads" \
        --outputPrefix="$prefix"

    # A wrapper can return successfully without producing a usable model.
    test -s "$prefix.rda"
    test -s "$prefix.varianceRatio.txt"
    tar -rf "$output" -C "$batch_tmp/models" \
        "$gene.rda" "$gene.varianceRatio.txt"
    rm -- "$prefix.rda" "$prefix.varianceRatio.txt" "$batch_tmp/$phenotype_file"
done
