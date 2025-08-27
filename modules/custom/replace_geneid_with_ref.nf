process REPLACE_GENEID_WITH_REF {
    tag "replace_geneid_with_ref"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d9bc3ff66b458498c00497ec45e762c23de475c90e6784fb393c16b472c77713/data' :
        'community.wave.seqera.io/library/python:3.10.18--6c172db5c7289169' }"

    input:
    path merged_gtf

    output:
    path "${merged_gtf.baseName}.refgene.gtf", emit: gtf
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    MERGED="$merged_gtf"
    # Derive the output name from the merged_gtf basename
    BASENAME=\$(basename "\$MERGED")
    OUT_NAME=\${BASENAME%.*}.refgene.gtf

    # Run the replacement script from the repository 'bin' folder
    replace_geneid_with_ref.py "\$MERGED" "\$OUT_NAME"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        replace_geneid_with_ref: \$(python3 --version 2>&1 | tr -d '\n')
    END_VERSIONS
    """
}
