process GFF_TO_GTF {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d9bc3ff66b458498c00497ec45e762c23de475c90e6784fb393c16b472c77713/data' :
        'community.wave.seqera.io/library/python:3.10.18--6c172db5c7289169' }"

    input:
    path plasmid_annotation
    val extension

    output:
    path "${plasmid_annotation.baseName}.fixed.gtf"

    script:
    """
    gff_to_gtf.py $plasmid_annotation -o ${plasmid_annotation.baseName}.fixed.gtf -e $extension
    """
}
