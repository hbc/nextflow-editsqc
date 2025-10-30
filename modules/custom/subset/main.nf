process SAMTOOLS_VIEW_REGION {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)
    path qname
    val index_format

    output:
    tuple val(meta), path("${prefix}.${meta.region}.${file_type}"),                                    emit: bam,              optional: true
    tuple val(meta), path("${prefix}.${meta.region}.${file_type}.bai"),                       emit: bai,              optional: true
    path  "versions.yml",                                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    file_type = args.contains("--output-fmt sam") ? "sam" :
                args.contains("--output-fmt bam") ? "bam" :
                args.contains("--output-fmt cram") ? "cram" :
                input.getExtension()

    output_file = index_format ? "${prefix}.${meta.region}.${file_type}##idx##${prefix}.${meta.region}.${file_type}.${index_format} --write-index" : "${prefix}.${meta.region}.${file_type}"
    // Can't choose index type of unselected file
    readnames = qname ? "--qname-file ${qname} --output-unselected ${prefix}.unselected.${file_type}": ""

    if (index_format) {
        if (!index_format.matches('bai|csi|crai')) {
            error "Index format not one of bai, csi, crai."
        } else if (file_type == "sam") {
            error "Indexing not compatible with SAM output"
        }
    }
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        view \\
        --threads ${task.cpus-1} \\
        ${reference} \\
        ${readnames} \\
        $args \\
        -o ${output_file} \\
        $input \\
        ${meta.region}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    file_type = args.contains("--output-fmt sam") ? "sam" :
                args.contains("--output-fmt bam") ? "bam" :
                args.contains("--output-fmt cram") ? "cram" :
                input.getExtension()
    default_index_format =
        file_type == "bam" ? "csi" :
        file_type == "cram" ? "crai" : ""
    index =  index_format ? "touch ${prefix}.${file_type}.${index_format}" : args.contains("--write-index") ? "touch ${prefix}.${file_type}.${default_index_format}" : ""
    unselected = qname ? "touch ${prefix}.unselected.${file_type}" : ""
    // Can't choose index type of unselected file
    unselected_index = qname && (args.contains("--write-index") || index_format) ? "touch ${prefix}.unselected.${file_type}.${default_index_format}" : ""

    if ("$input" == "${prefix}.${file_type}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    if (index_format) {
        if (!index_format.matches('bai|csi|crai')) {
            error "Index format not one of bai, csi, crai."
        } else if (file_type == "sam") {
            error "Indexing not compatible with SAM output."
        }
    }
    """
    touch ${prefix}.${file_type}
    ${index}
    ${unselected}
    ${unselected_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
