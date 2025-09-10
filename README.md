# Nextflow Nanopore Pipeline

## Introduction

This pipeline provides a reproducible workflow for processing Oxford Nanopore sequencing data using Nextflow and nf-core modules. It automates key steps from raw read alignment to transcript quantification and report generation, leveraging Docker for containerized execution and ensuring portability across systems.

## Pipeline Summary

The pipeline performs the following steps:

1. **Input Preparation**
   Accepts raw FASTQ files, reference genome FASTA, and annotation files (GFF/GTF).

2. **Read Alignment**
   Aligns Nanopore reads to the reference genome using Minimap2.

3. **BAM Processing**
   Sorts and indexes BAM files with Samtools.

4. **Annotation Conversion**
   Converts GFF annotations to GTF format and subsets annotations by region if specified.

5. **Transcriptome Assembly & Quantification**
   Assembles transcripts and quantifies expression using StringTie and Salmon.

6. **Quality Control**
   Performs QC on mapped reads and generates summary statistics.

7. **Report Generation**
   Creates IGV-compatible reports and summary files for downstream analysis.

## Supplying Inputs

Inputs are provided via Nextflow parameters, either in the command line or in a configuration file:

- `--fastq_files` : Path pattern to input FASTQ files (e.g., `data/*.fastq`)
- `--genome_fasta` : Path to reference genome FASTA file
- `--gene_annotation` : Path to genome annotation file (GFF or GTF)
- `--plasmid_fasta` : (Optional) Path to plasmid FASTA file
- `--plasmid_annotation` : (Optional) Path to plasmid annotation file
- `--region` : (Optional) Chromosome or region to subset annotation

Example command:
```bash
nextflow run main.nf --fastq_files 'data/*.fastq' --genome_fasta 'ref/genome.fa' --gene_annotation 'ref/genes.gff'
```

## Pipeline Outputs

The pipeline produces:

- **Aligned BAM files**: Sorted and indexed BAMs for each sample
- **Transcript GTF files**: Assembled transcript annotations
- **Quantification tables**: Gene and transcript expression estimates
- **Quality control reports**: Mapping statistics and QC summaries
- **IGV reports**: BED and annotation files for visualization in IGV
- **Log files**: Execution logs and pipeline metadata

All outputs are organized in the `results/` directory by default.

---

For more details on parameters and customization, see the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html) and the pipeline source code.
