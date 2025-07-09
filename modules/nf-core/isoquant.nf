#!/usr/bin/env nextflow

process ISOQUANT {
    container 'quay.io/biocontainers/isoquant:3.7.0--hdfd78af_0'

    input:
    tuple val(meta), path(fastq_files)
    path reference_genome
    path gene_annotation

    output:
    path 'isoquant_output_folder'

    script:
    """
    isoquant.py \
            --no_model_construction --gene_quantification all\
            --data_type nanopore \
            --reference $reference_genome \
            --genedb $gene_annotation \
            --fastq $fastq_files \
            --data_type nanopore -o isoquant_output_folder
    """
}
