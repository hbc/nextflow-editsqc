#!/usr/bin/env nextflow

nextflow.enable.dsl=2


include { ISOQUANT } from './modules/nf-core/isoquant'
include { TOULLIGQC } from './modules/nf-core/toulligqc'
include { TOULLIGQC as TOULLIGQC_GENOME } from './modules/nf-core/toulligqc'
include { GFFREAD as GFFREAD_FASTA} from './modules/nf-core/gffread'
include { GFFREAD as GFFREAD_GTF } from './modules/nf-core/gffread'
include { MINIMAP2_ALIGN } from './modules/nf-core/minimap2/align'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_GENOME } from './modules/nf-core/minimap2/align'
// include { SALMON_QUANT } from './modules/nf-core/salmon/quant'
// include { SALMON_QUANT as SALMON_QUANT_TX } from './modules/nf-core/salmon/quant'
// include { SALMON_INDEX } from './modules/nf-core/salmon/index'
include { STRINGTIE_STRINGTIE } from './modules/nf-core/stringtie/stringtie'
include { STRINGTIE_MERGE } from './modules/nf-core/stringtie/merge'
include { STRINGTIE_STRINGTIE as STRINGTIE_QUANT} from './modules/nf-core/stringtie/stringtie'
include { REPLACE_GENEID_WITH_REF } from './modules/custom/replace_geneid_with_ref'
include { GFF_TO_GTF as GFF_TO_GTF_PLASMID} from './modules/custom/gff_to_gtf'
include { GFF_TO_GTF as GFF_TO_GTF_GENOME} from './modules/custom/gff_to_gtf'
include { COMBINE_FASTA_ANNOTATION } from './modules/custom/combine_fasta_annotation'
include { FILTER_LONG_READS } from './modules/custom/filter_long_reads'
include { SAMTOOLS_VIEW }   from './modules/nf-core/samtools/view'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_REGION }   from './modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX }   from './modules/nf-core/samtools/index'
include { SAMTOOLS_IDXSTATS } from './modules/nf-core/samtools/idxstats'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TX }   from './modules/nf-core/samtools/index'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_TX } from './modules/nf-core/samtools/idxstats'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_TX } from './modules/nf-core/samtools/view'
include { PIPELINE_INFO } from './modules/custom/pipeline_info'

workflow {
    println(params)
    Channel.fromPath(params.genome_fasta, checkIfExists: true)
        .set { genome_fasta }
    Channel.fromPath(params.gene_annotation, checkIfExists: true)
        .set { gene_annotation }
    plasmid_annotation_ch = params.plasmid_annotation ? Channel.fromPath(params.plasmid_annotation, checkIfExists: true) : Channel.empty()
    plasmid_fasta_ch = params.plasmid_fasta ? Channel.fromPath(params.plasmid_fasta, checkIfExists: true) : Channel.empty()
    if (params.plasmid_annotation){
        plasmid_gtf = GFF_TO_GTF_PLASMID(plasmid_annotation_ch, params.extension)
    } else {
        plasmid_gtf = Channel.empty()
    }

    Channel.fromFilePairs(params.fastq_files, flat: true)
        .map { id, file -> tuple([id: id, name: file.baseName, single_end: params.single_end], file) }
        //.view()
        .set { fastq_files }

    gene_gtf = GFF_TO_GTF_GENOME(gene_annotation, params.extension)
    combined = COMBINE_FASTA_ANNOTATION(genome_fasta, plasmid_fasta_ch.ifEmpty([]), gene_gtf, plasmid_gtf.ifEmpty([]))
    // combined.out.view()
    gff_with_meta = combined.annotation
                        .map { file -> tuple([id: file.baseName, name: file.baseName], file) }
    fasta_with_meta = combined.fasta
                        .map { file -> tuple([id: file.baseName, name: file.baseName], file) }

    // get gtf with transcripts and exons
    gtf_with_txs = GFFREAD_GTF(gff_with_meta, combined.fasta)
    // gtf_with_txs.gtf
    //     .view()
    // get the transcript fasta file
    tx_output = GFFREAD_FASTA(gtf_with_txs.gtf, combined.fasta)
    // tx_output.gffread_fasta
    //     .view()

    tx_output.gffread_fasta
        .map { tuple -> tuple[1] }
        // .view()
        .set { fasta_alone }
    gtf_with_txs.gtf
        .map { tuple -> tuple[1] }
        // .view()
        .set { gff_alone }

    filtered_reads = FILTER_LONG_READS(fastq_files)

    TOULLIGQC(filtered_reads)

    // ISOQUANT(filtered_reads, combined.fasta, combined.annotation)

    align_output = MINIMAP2_ALIGN(filtered_reads, tx_output.gffread_fasta, "bam", "bai", false, true)
    // align_output.bam
    //     .view()

    // Map reads to genome fasta
    genome_align_output = MINIMAP2_ALIGN_GENOME(filtered_reads, fasta_with_meta, "bam", "bai", false, true)
    genome_idx = SAMTOOLS_INDEX(genome_align_output.bam)
    genome_bam = genome_align_output.bam
        .join(genome_idx.bai)
    SAMTOOLS_IDXSTATS(genome_bam)
    genome_bam
        // .view()

    if (params.region) {
        SAMTOOLS_VIEW_REGION(genome_bam, fasta_with_meta, [], 'bai')
    }
    only_mapped = SAMTOOLS_VIEW(genome_bam, fasta_with_meta, [], 'bai')
    TOULLIGQC_GENOME(only_mapped.bam)
    // Old Salmon quantification steps
    // Prepare gentrome and decoys for Salmon
    // salmon_index = SALMON_INDEX(combined.fasta, fasta_alone)
    // SALMON_QUANT(filtered_reads, salmon_index.index, combined.annotation, fasta_alone, false, 'A')
    //     tx_idx = SAMTOOLS_INDEX_TX(align_output.bam)
    // tx_bam = align_output.bam
    //     .join(tx_idx.bai)
    // SAMTOOLS_IDXSTATS_TX(tx_bam)
    // only_tx_mapped = SAMTOOLS_VIEW_TX(tx_bam, tx_output.gffread_fasta, [], 'bai')
    // SALMON_QUANT(align_output.bam, combined.annotation, fasta_alone, true, 'A')

    // Quantifie with stringtie
    stringtie_res = STRINGTIE_STRINGTIE(genome_align_output.bam, gff_alone)
    stringtie_res.transcript_gtf
        .map { tuple -> tuple[1] }
        // .view()
        .set { string_gtf_alone }

    stringtie_new = STRINGTIE_MERGE(string_gtf_alone, gff_alone)
    stringtie_fixed = REPLACE_GENEID_WITH_REF(stringtie_new.gtf)
    STRINGTIE_QUANT(genome_align_output.bam, stringtie_fixed.gtf)
    // Create a pipeline info file and publish it
    pipeline_info = PIPELINE_INFO()

    // Ensure it is available in the final workflow outputs (join to a channel used later)
    // pipeline_info
    //     .view()
    //     .set { final_pipeline_info }
}
