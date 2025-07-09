process GFF_TO_GTF {
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
