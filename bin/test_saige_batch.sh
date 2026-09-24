#!/usr/bin/env bash
# Run SAIGE Steps 2 and 3 sequentially for each gene; retain two batch tables.
set -euo pipefail
models=$1
vcf=$2
phenotype=$3
output_dir=$4
min_mac=$5
max_missing=$6
markers_per_chunk=$7

batch_tmp=$(mktemp -d "${PWD}/saige-test.XXXXXX")
trap 'rm -rf -- "$batch_tmp"' EXIT
tar -xf "$models" -C "$batch_tmp"
mkdir -p "$output_dir"
summary="$output_dir/gene_results.tsv"
printf 'phenotype_name\tgene_id\tstatus\tn_variants\tACAT_p\ttop_MarkerID\ttop_pval\n' > "$summary"

tail -n +2 "$batch_tmp/manifest.tsv" |
while IFS=$'\t' read -r gene chrom cis_start cis_end phenotype_file; do
    printf 'Testing %s\n' "$gene"
    # Reused for every gene and removed with the task's temporary directory.
    printf '%s\t%s\t%s\n' "$chrom" "$cis_start" "$cis_end" > "$batch_tmp/cis_region.tsv"
    associations="$batch_tmp/$gene.cis.tsv"
    step2_tests_qtl.R \
        --vcfFile="$vcf" --vcfField=GT \
        --GMMATmodelFile="$batch_tmp/$gene.rda" \
        --varianceRatioFile="$batch_tmp/$gene.varianceRatio.txt" \
        --rangestoIncludeFile="$batch_tmp/cis_region.tsv" --chrom="$chrom" \
        --minMAC="$min_mac" --minMAF=0 --maxMissing="$max_missing" \
        --is_imputed_data=FALSE --LOCO=FALSE --SPAcutoff=2 \
        --markers_per_chunk="$markers_per_chunk" --is_overwrite_output=TRUE \
        --SAIGEOutputFile="$associations"

    test -s "$associations"
    n_variants=$(awk 'NR > 1 && NF { n++ } END { print n+0 }' "$associations")
    if [ "$n_variants" -eq 0 ]; then
        printf '%s\t%s\tno_testable_variants\t0\tNA\tNA\tNA\n' \
            "$phenotype" "$gene" >> "$summary"
    else
        gene_pvalue="$batch_tmp/$gene.gene_pvalue.tsv"
        step3_gene_pvalue_qtl.R \
            --assocFile="$associations" --geneName="$gene" \
            --genePval_outputFile="$gene_pvalue"
        test "$(awk 'NR > 1 && NF { n++ } END { print n+0 }' "$gene_pvalue")" -eq 1
        awk -v phenotype="$phenotype" -v gene="$gene" -v n="$n_variants" \
            'BEGIN { FS = OFS = "\t" }
             NR > 1 && NF {
                 if ($1 != gene) exit 1
                 print phenotype, $1, "tested", n, $2, $3, $4
             }' "$gene_pvalue" >> "$summary"
        rm -- "$gene_pvalue"
    fi

    header=$(head -n 1 "$associations")
    if [ ! -e "$batch_tmp/associations.tsv" ]; then
        expected_header=$header
        printf 'phenotype_name\tgene_id\t%s\n' "$header" > "$batch_tmp/associations.tsv"
    elif [ "$header" != "$expected_header" ]; then
        printf 'Association columns changed for %s\n' "$gene" >&2
        exit 1
    fi
    awk -v phenotype="$phenotype" -v gene="$gene" \
        'BEGIN { OFS = "\t" } NR > 1 && NF { print phenotype, gene, $0 }' \
        "$associations" >> "$batch_tmp/associations.tsv"
    rm -- "$associations" "$batch_tmp/$gene.rda" "$batch_tmp/$gene.varianceRatio.txt"
done

gzip -c "$batch_tmp/associations.tsv" > "$output_dir/associations.tsv.gz"
