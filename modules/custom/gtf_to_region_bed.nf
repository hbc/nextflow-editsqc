process GTF_TO_REGION_BED {
    label 'process_low'

    input:
    path gtf
    val region

    output:
    path('regions.bed'), emit: bed

    script:
    """
    set -euo pipefail

   # If region is provided, subset GTF to that chromosome; else use all
    if [[ -n "${region}" && "${region}" != "null" ]]; then
        awk -v chr="${region}" '\$1 == chr' ${gtf} > subset.gtf
    else
        cp ${gtf} subset.gtf
    fi

    # Convert subset GTF to BED format
    awk '\$3 == "gene"' subset.gtf | \\
        awk 'BEGIN {OFS="\\t"} {print \$1, \$4-1, \$5, \$10}' | \\
        sed 's/"//g; s/;//g' > regions.bed

    echo "Generated regions.bed from subset.gtf"
    """
}
