process COMBINE_FASTA_ANNOTATION {
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/eb11a25875463a2bb2cca2dbef9103faabb1be227465732fa18e7373ee6eb66f/data' :
        'community.wave.seqera.io/library/seqtk:r93--b54ec2a2e8839010' }"

    input:
    path genome_fasta
    path plasmid_fasta
    path gene_annotation
    path plasmid_annotation

    output:
    path 'combined.fasta', emit: fasta
    path 'combined.gtf', emit: annotation

    script:
    """
    set -euo pipefail

    # Handle genome fasta decompression if needed
    GENOME_FILE="${genome_fasta}"
    if [[ "${genome_fasta}" == *.gz ]]; then
        echo "Decompressing genome fasta: ${genome_fasta}"
        gzip -cd "${genome_fasta}" > genome.uncompressed.fasta
        GENOME_FILE=genome.uncompressed.fasta
    fi

    # Handle plasmid fasta decompression if needed
    PLASMID_FILE=""
    if [[ "${plasmid_fasta}" != "NO_FILE" && -n "${plasmid_fasta}" ]]; then
        if [[ "${plasmid_fasta}" == *.gz ]]; then
            echo "Decompressing plasmid fasta: ${plasmid_fasta}"
            gzip -cd "${plasmid_fasta}" > plasmid.uncompressed.fasta
            PLASMID_FILE=plasmid.uncompressed.fasta
        else
            PLASMID_FILE="${plasmid_fasta}"
        fi
    fi

    # Concatenate FASTA files
    if [[ -n "\$PLASMID_FILE" ]]; then
        echo "Combining genome and plasmid FASTA files"
        cat "\$GENOME_FILE" "\$PLASMID_FILE" > tmp.fasta
    else
        echo "Using genome FASTA only"
        cat "\$GENOME_FILE" > tmp.fasta
    fi

    # Combine annotations
    if [[ "${plasmid_annotation}" != "NO_FILE" && -n "${plasmid_annotation}" ]]; then
        echo "Combining genome and plasmid annotations"
        cat "${gene_annotation}" "${plasmid_annotation}" > combined.gtf
    else
        echo "Using genome annotation only"
        cp "${gene_annotation}" combined.gtf
    fi

    # Format final FASTA with consistent line length
    seqtk seq -l 80 tmp.fasta > combined.fasta
    
    echo "Combined FASTA and GTF files created successfully"
    """
}