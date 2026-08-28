process PIPELINE_INFO {
    tag 'pipeline_info'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d9bc3ff66b458498c00497ec45e762c23de475c90e6784fb393c16b472c77713/data' :
        'community.wave.seqera.io/library/python:3.10.18--6c172db5c7289169' }"

    output:
    path 'methods.txt', emit: pipeline_info

    script:
    """
    set -euo pipefail

    cat > methods.txt <<'EOF'
    Pipeline: nextflow-editsqc
    Description: Workflow information and modes used for this run

    Tools and modes:
    - minimap2: unspliced
    - stringtie: reference-only
    - extension: ${params.extension}

    Generated: \$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    EOF
    """
}
