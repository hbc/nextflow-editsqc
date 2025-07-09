#!/usr/bin/env nextflow

nextflow.enable.dsl=2


include { TOULLIGQC } from './modules/nf-core/toulligqc'
include { ISOQUANT } from './modules/nf-core/isoquant'
include { GFFREAD } from './modules/nf-core/gffread'
include { MINIMAP2_ALIGN } from './modules/nf-core/minimap2/align'
include { SALMON_QUANT } from './modules/nf-core/salmon/quant'

workflow {
    Channel.fromPath(params.genome_fasta, checkIfExists: true)
        .set { genome_fasta }

    Channel.fromPath(params.gene_annotation, checkIfExists: true)
        .set { gene_annotation }

    Channel.fromFilePairs(params.fastq_files, flat: true)
        .map { id, file -> tuple([id: id, name: file.baseName], file) }
        //.view()
        .set { fastq_files }

    TOULLIGQC(fastq_files)

    ISOQUANT(fastq_files, genome_fasta, gene_annotation)

    gff_with_meta = gene_annotation
                        .map { file -> tuple([id: file.baseName, name: file.baseName], file) }

    gff_output = GFFREAD(gff_with_meta, genome_fasta)

    // gff_output.gffread_fasta
    //                     .map { file -> tuple([id: file.baseName, name: file.baseName], file) }
    //                     .view()
    //                     .set {fasta_with_meta}

    align_output = MINIMAP2_ALIGN(fastq_files, gff_output.gffread_fasta, true, '', false, true)

    gff_output.gffread_fasta
        .map { file -> tuple([id: file.baseName, name: file.baseName], file) }
        .set { fasta_alone }

    gff_output.gffread_fasta
        .map { fasta -> fasta[1] }
        //.view()
        .set { fasta_alone }

    SALMON_QUANT(align_output.bam, gene_annotation, fasta_alone, true, 'A')
}
