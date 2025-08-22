process COMBINE_FASTA_ANNOTATION {
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

    # Use Nextflow-injected path variables for inputs
    GENOME="$genome_fasta"
    PLASMID="$plasmid_fasta"

    # Decompress if gzipped
    if [[ "\$GENOME" == *.gz ]]; then
        echo "Decompressing genome fasta: \$GENOME"
        gzip -cd "\$GENOME" > genome.uncompressed.fasta
        GENOME=genome.uncompressed.fasta
    fi

    if [[ "\$PLASMID" == *.gz ]]; then
        echo "Decompressing plasmid fasta: \$PLASMID"
        gzip -cd "\$PLASMID" > plasmid.uncompressed.fasta
        PLASMID=plasmid.uncompressed.fasta
    fi

    # Concatenate and reformat
    cat "\$GENOME" "\$PLASMID" > tmp.fasta
    cat $gene_annotation $plasmid_annotation > combined.gtf
    seqtk seq -l 80 tmp.fasta > combined.fasta
    """
}
