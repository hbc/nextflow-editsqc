#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { MAIN_ANALYSIS } from './subworkflows/local/main_analysis'

workflow {
    // Create input channels
    // Channel.fromPath(params.genome_fasta, checkIfExists: true)
    //     .set { genome_fasta }

    // Channel.fromPath(params.gene_annotation, checkIfExists: true)
    //     .set { gene_annotation }

    // plasmid_annotation_ch = params.plasmid_annotation ?
    //     Channel.fromPath(params.plasmid_annotation, checkIfExists: true) :
    //     Channel.empty()

    // plasmid_fasta_ch = params.plasmid_fasta ?
    //     Channel.fromPath(params.plasmid_fasta, checkIfExists: true) :
    //     Channel.empty()

    // Channel.fromFilePairs(params.fastq_files, flat: true)
    //     .map { id, file -> tuple([id: id, name: file.baseName, single_end: params.single_end], file) }
    //     .set { fastq_files }

    // // Combine all inputs into a single tuple channel
    // input_tuple_ch = genome_fasta
    //     .combine(gene_annotation)
    //     .combine(plasmid_annotation_ch.ifEmpty([]))
    //     .combine(plasmid_fasta_ch.ifEmpty([]))
    //     .combine(fastq_files)
    //     //.combine(Channel.value(params.valid_features))
    //     //.combine(Channel.value(params.extension))
    //     //.combine(Channel.value(params.single_end))
    //     .combine(Channel.value(params.region))
    //     .view()

// Read the samplesheet and create input channels
    Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            // Create metadata map
            def meta = [id: row.sample, name: row.sample, single_end: params.single_end, region: row.region]

            // Process fastq files - handle both single files and file pairs
            def fastq_file = file(row.fastq_files, checkIfExists: true)
            def fastq_tuple = tuple(meta, fastq_file)

            // Create file objects for required files
            def genome_fasta = file(row.genome_fasta, checkIfExists: true)
            def gene_annotation = file(row.gene_annotation, checkIfExists: true)

            // Handle optional plasmid files
            def plasmid_fasta = row.plasmid_fasta && row.plasmid_fasta != '' ?
                file(row.plasmid_fasta, checkIfExists: true) : null
            def plasmid_annotation = row.plasmid_annotation && row.plasmid_annotation != '' ?
                file(row.plasmid_annotation, checkIfExists: true) : null

            // Create the input tuple matching your current structure
            return tuple(
                genome_fasta,
                gene_annotation,
                plasmid_annotation,
                plasmid_fasta,
                fastq_tuple,
                row.region
            )
        }
        .set { input_tuple_ch }

    // View the processed input tuples
    // input_tuple_ch.view()
    // Run the main analysis subworkflow with single tuple input
    MAIN_ANALYSIS(input_tuple_ch, params.valid_features, params.extension, params.single_end)

    // Access outputs from the subworkflow
    MAIN_ANALYSIS.out.combined_fasta.view { "Combined FASTA: $it" }
    MAIN_ANALYSIS.out.pipeline_info.view { "Pipeline info: $it" }
}
