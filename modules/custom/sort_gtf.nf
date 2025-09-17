process SORT_GTF {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d9bc3ff66b458498c00497ec45e762c23de475c90e6784fb393c16b472c77713/data' :
        'community.wave.seqera.io/library/python:3.10.18--6c172db5c7289169' }"

    input:
    path gtf

    output:
    path "*_sorted.gtf", emit: gtf

    script:
    """
    set -euo pipefail   

    sort -k1,1 -k4,4n $gtf > ${gtf.baseName}_sorted.gtf
    """
}