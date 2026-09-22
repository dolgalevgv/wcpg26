#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process PREPARE_DONORS {
    tag "${donors.baseName}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd4d58047a5a9c36d1ab9f6b9ac3242d1ee123e773f96548a3ddf5a26fb15453/data'

    input:
    path(donors)

    output:
    path("${donors.baseName}_vcf.txt"), emit: donors_vcf

    script:
    """
    prepare_donors.py ${donors} ${donors.baseName}_vcf.txt
    """
}

process BCFTOOLS_EXTRACT_DONORS {
    tag "${vcf.getBaseName(2)}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a08933da914fc6b3650dbe842d7c2e4d95155d9a57a8c49c3756ad7248b7fc2c/data'

    input:
    path(vcf)
    path(donors_vcf)
    
    output:
    path("${vcf.getBaseName(2)}_donors.vcf.gz")

    script:
    """
    bcftools view ${vcf} \\
        --samples-file <(cut -d' ' -f1 ${donors_vcf}) \\
        --output-type z \\
        --output ${vcf.getBaseName(2)}_extracted.vcf.gz
    
    bcftools reheader ${vcf.getBaseName(2)}_extracted.vcf.gz \\
        --samples ${donors_vcf} \\
        --output ${vcf.getBaseName(2)}_donors.vcf.gz
    """
}

process BCFTOOLS_LIFTOVER {
    tag "${vcf.getBaseName(2)}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a08933da914fc6b3650dbe842d7c2e4d95155d9a57a8c49c3756ad7248b7fc2c/data'

    input:
    path(vcf)
    path(source_fasta)
    path(target_fasta)
    path(chain_file)

    output:
    path("${vcf.getBaseName(2)}_lift_sorted.vcf.gz"), emit: lifted
    path("${vcf.getBaseName(2)}_lift_rejected_sorted.vcf.gz"), emit: rejected

    script:
    """
    bcftools norm ${vcf} \\
        --check-ref s \\
        --do-not-normalize \\
        --fasta-ref ${source_fasta} \\
        --output-type u \\
        | bcftools +liftover \\
            --output-type u \\
            -- \\
            --src-fasta-ref ${source_fasta} \\
            --fasta-ref ${target_fasta} \\
            --chain ${chain_file} \\
            --write-reject \\
            --reject ${vcf.getBaseName(2)}_lift_rejected_unsorted.vcf.gz \\
            --reject-type z \\
        | bcftools sort \\
            --output-type z \\
            --output ${vcf.getBaseName(2)}_lift_sorted.vcf.gz

    bcftools sort ${vcf.getBaseName(2)}_lift_rejected_unsorted.vcf.gz \\
        --output-type z \\
        --output ${vcf.getBaseName(2)}_lift_rejected_sorted.vcf.gz
    """
}

process BCFTOOLS_VARIANT_QC {
    tag "${vcf.getBaseName(2)}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a08933da914fc6b3650dbe842d7c2e4d95155d9a57a8c49c3756ad7248b7fc2c/data'

    input:
    path(vcf)

    output:
    path("${vcf.getBaseName(2)}_filt.vcf.gz"), emit: vcf

    script:
    """
    bcftools view ${vcf} \\
        --min-ac 10:minor \\
        --include 'F_MISSING <= 0.02' \\
        --min-alleles 2 \\
        --max-alleles 2 \\
        --types snps \\
        --output-type z \\
        --output ${vcf.getBaseName(2)}_filt.vcf.gz
    """
}

process BCFTOOLS_INDEX {
    tag "${vcf.getBaseName(2)}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a08933da914fc6b3650dbe842d7c2e4d95155d9a57a8c49c3756ad7248b7fc2c/data'

    input:
    path(vcf)

    output:
    tuple path("${vcf}"), path("${vcf}.csi"), emit: vcf

    script:
    """
    bcftools index --csi ${vcf}
    """
}

process PLINK_MAKE_PFILE {
    tag "${vcf.getBaseName(2)}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7fbbbd635adc17f214e69145009a0d1d0411c350b5e70eb14b5aa68d79a3fa1b/data'

    input:
    path(vcf)

    output:
    tuple path("${vcf.getBaseName(2)}.pgen"), path("${vcf.getBaseName(2)}.psam"), path("${vcf.getBaseName(2)}.pvar"), emit: pfile

    script:
    """
    plink2 \\
        --vcf ${vcf} \\
        --set-all-var-ids '@:#:\$r:\$a' \\
        --make-pgen \\
        --out ${vcf.getBaseName(2)}
    """
}

process PLINK_INDEP_PAIRWISE {
    tag "${pgen.baseName}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7fbbbd635adc17f214e69145009a0d1d0411c350b5e70eb14b5aa68d79a3fa1b/data'

    input:
    tuple path(pgen), path(psam), path(pvar)
    path(ld_exclude_bed)

    output:
    tuple path(pgen), path(psam), path(pvar), path("${pgen.baseName}.prune.in"), emit: pfile

    script:
    """
    plink2 \\
        --pfile ${pgen.baseName} \\
        --exclude bed1 ${ld_exclude_bed} \\
        --bad-ld \\
        --indep-pairwise 200kb 0.2 \\
        --out ${pgen.baseName}
    """
}

process PLINK_PCA {
    tag "${pgen.baseName}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7fbbbd635adc17f214e69145009a0d1d0411c350b5e70eb14b5aa68d79a3fa1b/data'

    publishDir "${params.outdir}/plink2", copy: true, pattern: "*.eigenvec"
    input:
    tuple path(pgen), path(psam), path(pvar), path(pruned)

    output:
    path("${pgen.baseName}.eigenvec")

    script:
    """
    plink2 \
        --pfile ${pgen.baseName} \\
        --extract ${pruned} \\
        --bad-freqs \\
        --pca ${params.n_pcs} \\
        --out ${pgen.baseName}
    """
}

process PREPARE_GENE_REGIONS {
    tag "${gtf.baseName}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd4d58047a5a9c36d1ab9f6b9ac3242d1ee123e773f96548a3ddf5a26fb15453/data'
    publishDir "${params.outdir}/gene_regions", mode: 'copy'
    
    input:
    path(gtf)
    path(fai)

    output:
    path("gene_regions.tsv"), emit: regions

    script:
    """
    prepare_gene_regions.py ${gtf} ${fai} gene_regions.tsv ${params.cis_window}
    """
}

process PREPARE_QTL_PHENOTYPES {
    tag "${adata.baseName}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5c688ea7f743de8aa32394006d13bd87e1d5d03df3e1e3b0443907de1ac786c9/data'

    publishDir "${params.outdir}/qtl_phenotypes/", mode: 'copy', pattern: "{*.csv,**/*.csv}"

    input:
    path(adata)
    path(donors)
    path(pca)
    path(regions)

    output:
    path("manifest.csv"), emit: manifest
    path("*_*"), type: 'dir', optional: true, emit: phenotypes

    script:
    """
    prepare_qtl_phenotypes.py \\
        "${adata}" \\
        "${donors}" \\
        "${pca}" \\
        "${regions}" \\
        --cell_type_col ${params.cell_type_col} \\
        --treat_col ${params.treat_col} \\
        --qcovar ${params.qcovar} \\
        --min-n-cells ${params.min_n_cells} \\
        --min-n-donors ${params.min_n_donors} \\
        --min-cells-frac ${params.min_cells_frac} \\
        --min-donor-frac ${params.min_donor_frac}
    """
}

process PLINK_PREPARE_SAIGE_BFILE {
    tag "${phenotype_name}"
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7f/7fbbbd635adc17f214e69145009a0d1d0411c350b5e70eb14b5aa68d79a3fa1b/data'

    publishDir { "${params.outdir}/saige_genotypes/${phenotype_name}" }, mode: 'copy', pattern: "*saige_vr*"

    input:
    tuple val(phenotype_name), path(phenotype_dir)
    tuple path(pgen), path(psam), path(pvar), path(pruned)

    output:
    tuple val(phenotype_name), 
          path(phenotype_dir), 
          path("${phenotype_name}_saige_vr.bed"), 
          path("${phenotype_name}_saige_vr.fam"), 
          path("${phenotype_name}_saige_vr.bim"), 
          emit: qtl_inputs

    script:
    """
    tail -n +2 ${phenotype_dir}/donors.csv | cut -d, -f1 > donors_keep.txt

    plink2 \\
        --pfile ${pgen.baseName} \\
        --keep donors_keep.txt \\
        --extract ${pruned} \\
        --mac 20 \\
        --make-bed \\
        --out ${phenotype_name}_saige_vr
    """
}

workflow {
    vcf_ch = Channel.fromPath(params.vcf, checkIfExists: true)
    donors = file(params.donors, checkIfExists: true)
    source_fasta = file(params.source_fasta, checkIfExists: true)
    target_fasta = file(params.target_fasta, checkIfExists: true)
    target_fai = file("${params.target_fasta}.fai", checkIfExists: true)
    chain_file = file(params.chain_file, checkIfExists: true)
    ld_exclude_bed = file(params.ld_exclude_bed, checkIfExists: true)
    gene_gtf = file(params.gene_gtf, checkIfExists: true)
    adata = file(params.adata, checkIfExists: true)

    PREPARE_DONORS(donors)
    BCFTOOLS_EXTRACT_DONORS(vcf_ch, PREPARE_DONORS.out.donors_vcf)
    BCFTOOLS_LIFTOVER(
        BCFTOOLS_EXTRACT_DONORS.out,
        source_fasta,
        target_fasta,
        chain_file
    )
    BCFTOOLS_VARIANT_QC(BCFTOOLS_LIFTOVER.out.lifted)

    BCFTOOLS_INDEX(BCFTOOLS_VARIANT_QC.out.vcf)
    PLINK_MAKE_PFILE(BCFTOOLS_VARIANT_QC.out.vcf)
    PLINK_INDEP_PAIRWISE(PLINK_MAKE_PFILE.out.pfile, ld_exclude_bed)
    PLINK_PCA(PLINK_INDEP_PAIRWISE.out.pfile)

    PREPARE_GENE_REGIONS(gene_gtf, target_fai)

    PREPARE_QTL_PHENOTYPES(
        adata,
        donors,
        PLINK_PCA.out,
        PREPARE_GENE_REGIONS.out.regions
    )

    qtl_phenotypes_ch = PREPARE_QTL_PHENOTYPES.out.phenotypes
        .flatten()
        .map { phenotypes_dir ->
            tuple(phenotypes_dir.name, phenotypes_dir)
        }

    PLINK_PREPARE_SAIGE_BFILE(
        qtl_phenotypes_ch, 
        PLINK_INDEP_PAIRWISE.out.pfile.first()
        )
}
