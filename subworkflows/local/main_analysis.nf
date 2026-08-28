#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TOULLIGQC } from '../../modules/nf-core/toulligqc'
include { TOULLIGQC as TOULLIGQC_GENOME } from '../../modules/nf-core/toulligqc'
include { GFFREAD as GFFREAD_FASTA} from '../../modules/nf-core/gffread'
include { GFFREAD as GFFREAD_GTF } from '../../modules/nf-core/gffread'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_GENOME } from '../../modules/nf-core/minimap2/align'
include { STRINGTIE_STRINGTIE } from '../../modules/nf-core/stringtie/stringtie'
include { STRINGTIE_MERGE } from '../../modules/nf-core/stringtie/merge'
include { STRINGTIE_STRINGTIE as STRINGTIE_QUANT} from '../../modules/nf-core/stringtie/stringtie'
include { REPLACE_GENEID_WITH_REF } from '../../modules/custom/replace_geneid_with_ref'
include { GFF_TO_GTF as GFF_TO_GTF_PLASMID} from '../../modules/custom/gff_to_gtf'
include { GFF_TO_GTF as GFF_TO_GTF_GENOME} from '../../modules/custom/gff_to_gtf'
include { GTF_TO_REGION_BED } from '../../modules/custom/gtf_to_region_bed'
include { COMBINE_FASTA_ANNOTATION } from '../../modules/custom/combine_fasta_annotation'
include { FILTER_LONG_READS } from '../../modules/custom/filter_long_reads'
include { SAMTOOLS_VIEW }   from '../../modules/nf-core/samtools/view'
include { SAMTOOLS_VIEW_REGION }   from '../../modules/custom/subset'
include { SAMTOOLS_INDEX }   from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_FAIDX }   from '../../modules/nf-core/samtools/faidx'
include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/samtools/idxstats'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TX }   from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_IDXSTATS as SAMTOOLS_IDXSTATS_TX } from '../../modules/nf-core/samtools/idxstats'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_TX } from '../../modules/nf-core/samtools/view'
include { IGVREPORTS } from '../../modules/nf-core/igvreports/main'
include { SORT_GTF } from '../../modules/custom/sort_gtf'
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'
include { PIPELINE_INFO } from '../../modules/custom/pipeline_info'

workflow MAIN_ANALYSIS {
    take:
    input_tuple  // channel: tuple of [genome_fasta, gene_annotation, plasmid_annotation, plasmid_fasta, fastq_files, valid_features, extension, single_end, region]
    valid_features
    extension
    single_end

    main:
    // Extract individual components from the input tuple
    genome_fasta = input_tuple.map { it[0] }
    gene_annotation = input_tuple.map { it[1] }
    plasmid_annotation = input_tuple.map { it[2] }
    plasmid_fasta = input_tuple.map { it[3] }
    fastq_files = input_tuple.map { it[4] }
    // valid_features = input_tuple.map { it[5] }
    // extension = input_tuple.map { it[6] }
    // single_end = input_tuple.map { it[7] }
    region = input_tuple.map { it[5] }
    region = region.ifEmpty("[]")
    // Convert plasmid annotation if provided
    if (plasmid_annotation) {
        plasmid_gtf = GFF_TO_GTF_PLASMID(plasmid_annotation, valid_features, extension)
    } else {
        plasmid_gtf = Channel.empty()
    }

    // Convert genome annotation
    gene_gtf = GFF_TO_GTF_GENOME(gene_annotation, valid_features, extension)

    // Combine FASTA and annotation files
    combined = COMBINE_FASTA_ANNOTATION(
        genome_fasta,
        plasmid_fasta.ifEmpty([]),
        gene_gtf,
        plasmid_gtf.ifEmpty([])
    )

    // Sort, zip, and index combined GTF
    combined_sorted = SORT_GTF(combined.annotation)
    combined_sorted_with_meta = combined_sorted.gtf
                        .map { file -> tuple([id: file.baseName, name: file.baseName], file) }
    TABIX_BGZIPTABIX(combined_sorted_with_meta)

    // Prepare channels with metadata
    gff_with_meta = combined.annotation
                        .map { file -> tuple([id: file.baseName, name: file.baseName], file) }
    fasta_with_meta = combined.fasta
                        .map { file -> tuple([id: file.baseName, name: file.baseName], file) }

    // Index the combined FASTA
    SAMTOOLS_FAIDX(fasta_with_meta, [[],[]], false)

    // Get GTF with transcripts and exons
    gtf_with_txs = GFFREAD_GTF(gff_with_meta, combined.fasta)

    // Get the transcript FASTA file
    tx_output = GFFREAD_FASTA(gtf_with_txs.gtf, combined.fasta)

    // Extract FASTA and GTF without metadata
    fasta_alone = tx_output.gffread_fasta
        .map { tuple -> tuple[1] }

    gff_alone = gtf_with_txs.gtf
        .map { tuple -> tuple[1] }

    // Filter long reads
    filtered_reads = FILTER_LONG_READS(fastq_files)

    // Quality control of filtered reads
    TOULLIGQC(filtered_reads)

    // Align reads to transcriptome
    align_output = MINIMAP2_ALIGN(filtered_reads, tx_output.gffread_fasta, "bam", "bai", false, true)

    // Align reads to genome
    genome_align_output = MINIMAP2_ALIGN_GENOME(filtered_reads, fasta_with_meta, "bam", "bai", false, true)

    // Index genome alignments
    genome_idx = SAMTOOLS_INDEX(genome_align_output.bam)
    genome_bam = genome_align_output.bam
        .join(genome_idx.bai)

    // Get alignment statistics
    SAMTOOLS_IDXSTATS(genome_bam)

    // Optional region-specific analysis.
    // Samples with a blank `region` column are dropped here, so every downstream
    // region process simply gets no work rather than running on a missing contig.
    region_bam = genome_bam.filter { meta, bam, bai -> meta.region }

    // Per-sample reference files, keyed by sample id so they can be joined on meta
    // instead of cross-multiplied with .combine()
    ref_by_sample = fastq_files
        .map { meta, reads -> meta.id }
        .merge(combined.fasta)
        .merge(combined.annotation)
        .map { id, fasta, gtf -> tuple(id, fasta, gtf) }

    // View specific region
    view_region = SAMTOOLS_VIEW_REGION(region_bam, fasta_with_meta, [], 'bai')

    // Create BED file for the region of each sample
    region_bed = GTF_TO_REGION_BED(
        region_bam
            .map { meta, bam, bai -> tuple(meta.id, meta) }
            .join(ref_by_sample)
            .map { id, meta, fasta, gtf -> tuple(meta, gtf, meta.region) }
    )

    // Assemble IGV inputs by joining on meta, so each sample keeps its own
    // bed / bam / annotation instead of pairing with another sample's files
    igv_input_ch = region_bed.bed
        .join(view_region.bam)
        .join(view_region.bai)
        .map { meta, bed, bam, bai -> tuple(meta.id, meta, bed, bam, bai) }
        .join(ref_by_sample)
        .map { id, meta, bed, bam, bai, fasta, gtf ->
            tuple(meta, bed, [gtf, bam], [bai])
        }

    fasta_input_ch = region_bam
        .map { meta, bam, bai -> tuple(meta.id, meta) }
        .join(ref_by_sample)
        .map { id, meta, fasta, gtf -> tuple(meta, fasta, file("${fasta}.fai")) }

    // Generate IGV reports
    IGVREPORTS(igv_input_ch, fasta_input_ch)

    // Filter for mapped reads only
    only_mapped = SAMTOOLS_VIEW(genome_bam, fasta_with_meta, [], 'bai')

    // Quality control of mapped reads
    TOULLIGQC_GENOME(only_mapped.bam)

    // Quantification with StringTie
    stringtie_res = STRINGTIE_STRINGTIE(genome_align_output.bam, gff_alone)

    string_gtf_alone = stringtie_res.transcript_gtf
        .map { tuple -> tuple[1] }

    // Merge StringTie results
    stringtie_new = STRINGTIE_MERGE(string_gtf_alone, gff_alone)

    // Replace gene IDs with reference
    stringtie_fixed = REPLACE_GENEID_WITH_REF(stringtie_new.gtf)

    // Final quantification
    STRINGTIE_QUANT(genome_align_output.bam, stringtie_fixed.gtf)

    // Create pipeline info
    pipeline_info = PIPELINE_INFO()

    emit:
    combined_fasta           = combined.fasta
    combined_annotation      = combined.annotation
    filtered_reads          = filtered_reads
    genome_alignments       = genome_align_output.bam
    transcriptome_alignments = align_output.bam
    stringtie_results       = STRINGTIE_QUANT.out.abundance
    pipeline_info           = pipeline_info
    // toulligqc_reports       = TOULLIGQC.out
    // toulligqc_genome_reports = TOULLIGQC_GENOME.out
    // idxstats                = SAMTOOLS_IDXSTATS.out
}
