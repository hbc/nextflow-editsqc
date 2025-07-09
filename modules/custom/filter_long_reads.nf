process FILTER_LONG_READS {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/55/5517fc4f67c27658bbbff2b18da31b88a8572346e90075739450df40d6b5e315/data' :
        'community.wave.seqera.io/library/fastplong:0.3.0--7ac22bc1d1bce90c' }"

    input:
    tuple val(meta), path(reads)


    output:
    tuple val(meta), path("${reads.baseName}_filtered.fastq.gz"), emit: filtered_reads

    script:
    def args = task.ext.args ?: ''
    """
    fastplong filter -i $reads -o ${reads.baseName}_filtered.fastq.gz $args
    """
}
