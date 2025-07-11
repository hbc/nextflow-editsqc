#!/usr/bin/env nextflow

nextflow.enable.dsl=2


include { TOULLIGQC } from './modules/nf-core/toulligqc'
include { TOULLIGQC as TOULLIGQC_GENOME } from './modules/nf-core/toulligqc'
include { ISOQUANT } from './modules/nf-core/isoquant'
include { GFFREAD as GFFREAD_FASTA} from './modules/nf-core/gffread'
include { GFFREAD as GFFREAD_GTF } from './modules/nf-core/gffread'
include { MINIMAP2_ALIGN } from './modules/nf-core/minimap2/align'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_GENOME } from './modules/nf-core/minimap2/align'
include { SALMON_QUANT } from './modules/nf-core/salmon/quant'
include { GFF_TO_GTF as GFF_TO_GTF_PLASMID} from './modules/custom/gff_to_gtf'
include { GFF_TO_GTF as GFF_TO_GTF_GENOME} from './modules/custom/gff_to_gtf'
include { COMBINE_FASTA_ANNOTATION } from './modules/custom/combine_fasta_annotation'
include { FILTER_LONG_READS } from './modules/custom/filter_long_reads'
include { SAMTOOLS_VIEW }   from './modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX }   from './modules/nf-core/samtools/index'
include { SAMTOOLS_IDXSTATS } from './modules/nf-core/samtools/idxstats'

workflow {
    println(params)
    Channel.fromPath(params.genome_fasta, checkIfExists: true)
        .set { genome_fasta }
    Channel.fromPath(params.gene_annotation, checkIfExists: true)
        .set { gene_annotation }
    Channel.fromPath(params.plasmid_annotation, checkIfExists: true)
        .set { plasmid_annotation }
    Channel.fromPath(params.plasmid_fasta, checkIfExists: true)
        .set { plasmid_fasta }

    Channel.fromFilePairs(params.fastq_files, flat: true)
        .map { id, file -> tuple([id: id, name: file.baseName], file) }
        //.view()
        .set { fastq_files }

    plasmid_gtf = GFF_TO_GTF_PLASMID(plasmid_annotation, 50)
    gene_gtf = GFF_TO_GTF_GENOME(gene_annotation, 50)

    combined = COMBINE_FASTA_ANNOTATION(plasmid_fasta, genome_fasta, plasmid_gtf, gene_gtf)
    // combined.out.view()
    gff_with_meta = combined.annotation
                        .map { file -> tuple([id: file.baseName, name: file.baseName], file) }
    fasta_with_meta = combined.fasta
                        .map { file -> tuple([id: file.baseName, name: file.baseName], file) }
    // combined.fasta
    //     .view()
    // combined.annotation
    //     .view()
    // get gtf with transcripts and exons
    gtf_with_txs = GFFREAD_GTF(gff_with_meta, combined.fasta)
    // gtf_with_txs.gtf
    //     .view()
    // get the transcript fasta file
    tx_output = GFFREAD_FASTA(gtf_with_txs.gtf, combined.fasta)
    // tx_output.gffread_fasta
    //     .view()
    filtered_reads = FILTER_LONG_READS(fastq_files)

    TOULLIGQC(filtered_reads)

    ISOQUANT(filtered_reads, combined.fasta, combined.annotation)

    align_output = MINIMAP2_ALIGN(filtered_reads, tx_output.gffread_fasta, true, '', false, true)
    // align_output.bam
    //     .view()
    
    // Map reads to genome fasta
    genome_align_output = MINIMAP2_ALIGN_GENOME(filtered_reads, fasta_with_meta, true, '', false, true)
    genome_idx = SAMTOOLS_INDEX(genome_align_output.bam)
    genome_bam = genome_align_output.bam
        .join(genome_idx.bai)
    SAMTOOLS_IDXSTATS(genome_bam)
    genome_bam
        .view()
    only_mapped = SAMTOOLS_VIEW(genome_bam, fasta_with_meta, [], 'bai')
    TOULLIGQC_GENOME(only_mapped.bam)

    tx_output.gffread_fasta
        .map { tuple -> tuple[1] }
        // .view()
        .set { fasta_alone }

    SALMON_QUANT(align_output.bam, combined.annotation, fasta_alone, true, 'A')

}
