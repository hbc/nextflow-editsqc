process PIPELINE_INFO {
    tag 'pipeline_info'

    // No container needed; simple shell job

    output:
    path 'methods.txt', emit: pipeline_info

    script:
    """
    set -euo pipefail

    cat > methods.txt <<'EOF'
    Pipeline: nextflow-nanopore
    Description: Workflow information and modes used for this run

    Tools and modes:
    - minimap2: unspliced
    - stringtie: reference-only
    - extension: ${params.extension}

    Generated: \$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    EOF
    """
}
