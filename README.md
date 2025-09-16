# Nextflow Nanopore Pipeline

## Introduction

This pipeline provides a reproducible workflow for processing Oxford Nanopore sequencing data using Nextflow and nf-core modules. It automates key steps from raw read alignment to transcript quantification and report generation, leveraging Docker for containerized execution and ensuring portability across systems.

## Pipeline Summary

The pipeline performs the following steps:

1. **Input Preparation**
   Accepts raw FASTQ files, reference genome FASTA, and annotation files (GFF/GTF).

2. **Read Alignment**
   Aligns Nanopore reads to the reference genome using Minimap2.


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


### Specific details about FASTA and GTF

Below are specific formatting rules and examples to ensure your FASTA and GFF/GTF inputs are compatible with the pipeline.

- FASTA (genome or transcriptome)
  - Headers must start with a single '>' character followed by a unique sequence identifier (sequence IDs are used to match annotation seqids). Avoid spaces in the primary ID; use underscores if needed.
    - Good: `>chromosome` or `>plasmid_1` or `>TAC.gene01.1`
    - Avoid: `>chromosome 1 description` (the space makes the full header the ID)
  - The sequence identifier in the FASTA header must exactly match the seqid fields used in your GFF/GTF (case-sensitive).
  - Sequence lines may be wrapped; the pipeline will read standard FASTA wrapping.

- GFF3 (annotation)
  - The pipeline accepts GFF3. For best results follow these rules:
    - Seqids (first column) must match FASTA headers exactly (e.g., `chromosome`, `plasmid`, or transcript IDs used in transcriptome FASTA).
    - Use 1-based inclusive coordinates (standard GFF convention).
    - For GFF3: ensure feature lines include `Name`, `ID` and `Parent` attributes where appropriate (gene/mRNA/exon). Example GFF3 line:
      `plasmid\tmaker\tmRNA\t2501\t3500\t.\t+\t.\tID=TAC.plasmid_gene01.1;Name=gene1;Parent=geneX;`
  - If your GFF contains attributes or naming conventions that do not include `gene_id`/`transcript_id`, the pipeline's `gffread` and conversion modules attempt to translate common patterns, but explicit attributes are more reliable.

- Common pitfalls and tips
  - Mismatched seqids (FASTA header vs GFF seqid column) are the most frequent source of problems — check these first if features or mapping counts are missing.
  - Ensure no duplicate sequence IDs in FASTA; tools (StringTie) require unique IDs.
  - GFF3 multi-line attributes or unusual quoting can sometimes break parsers; provide standard, single-line attribute fields.
  - Provide a GFF3 — if you have GFF3 the pipeline will convert to GTF for StringTie where needed.

## Running the pipeline

Example command:
```bash
nextflow run main.nf --fastq_files 'data/*.fastq' --genome_fasta 'ref/genome.fa' --gene_annotation 'ref/genes.gff'
```

Or create a config file (`mysample.config`) and use it in the command line:

```
params{
  fastq_files= 'nanopore_np.fastq.gz'
  genome_fasta= './plasmids/plasmid1.fasta'
  gene_annotation= './plasmids/plasmid1.annotations.gff3'
  // plasmid_fasta= './genome_plasmid.fa'
  // plasmid_annotation= './genome_plasmid.gff'
  outdir= './results'
  extension = '0'
  region = "plasmid_chromosome_name"
}
```

and then run with:

```
nextflow run main.nf -profile docker -c mysample.config --outdir myresults
```

## Pipeline Outputs

The pipeline writes all results to the `results/` directory (or a custom `--outdir`). Below is a detailed description of the directories and files produced by this pipeline (examples taken from the included `results/` folder).

### Most useful files

- `combine_fasta_annotation/`
  - `combined.fasta` — combined transcriptome/genome FASTA used for transcript-level workflows (e.g. Salmon index/quantification).
  - `combined.gtf` — combined GTF annotation generated from input GFF/GTF files and any plasmid annotations.

- `minimap2/`
  - `genome/`
    - `nanopore_reads.bam` — reads aligned to the genome (sorted BAM).
    - `nanopore_reads.bam.bai` — BAM index for quick access.
    - `nanopore_reads.idxstats` — per-reference mapping counts (chromosome/plasmid read counts).
    - `qc/` — per-sample QC reports generated by ToulligQC (HTML report and images).
    - `region/` — optional subdirectories if reads were aligned or subset by genomic region.
  - `tx/`
    - `nanopore_reads.bam` — reads aligned to the transcriptome (if generated).
    - `nanopore_reads.bam.bai` — BAM index for transcriptome alignments.
    - `nanopore_reads.idxstats` — mapping counts per transcript.


- `stringtie_new/`
  - `nanopore_reads.transcripts.gtf` — transcripts assembled by StringTie from the sample BAM.
  - `nanopore_reads.coverage.gtf` — coverage tracks produced by StringTie (`-C` option) showing assembled loci coverage.
  - `nanopore_reads.gene.abundance.txt` — gene-level abundance table with Coverage, FPKM, and TPM columns.

### Complementary output files

- `filtered_reads/`
  - `nanopore_reads.fastq_filtered.fastq.gz` — input FASTQ after optional filtering (length/quality) performed by the `filter_long_reads` module.

- `gffread/`
  - `tx.gtf` — transcript GTF produced by converting and/or subsetting the original annotation (used by StringTie and Salmon).
  - `tx.seq.fasta` — transcriptome FASTA produced from the annotation and reference, used for transcript-level quantification.
  - `versions.yml` — versions of tools used to generate these files.

- `pipeline_info/`
  - `execution_report_<timestamp>.html` — full Nextflow execution report (runtime, resource usage, and process summaries).
  - `execution_timeline_<timestamp>.html` — timeline of workflow execution for debugging and review.
  - `execution_trace_<timestamp>.txt` — raw trace of tasks run by Nextflow.
  - `pipeline_dag_<timestamp>.html` — graphical representation of the pipeline DAG.

- `stringtie/` and `stringtie_new/`
  - `nanopore_reads.transcripts.gtf` — transcripts assembled by StringTie from the sample BAM.
  - `nanopore_reads.coverage.gtf` — coverage tracks produced by StringTie (`-C` option) showing assembled loci coverage.
  - `nanopore_reads.gene.abundance.txt` — gene-level abundance table with Coverage, FPKM, and TPM columns.
  - `stringtie.merged.gtf` / `stringtie.merged.refgene.gtf` — merged transcriptomes across samples and (optionally) reconciled with reference annotation.
  - `nanopore_reads.ballgown/` — Ballgown-compatible output for downstream expression analysis.

- `toulligqc/`
  - `nanopore_readsToulligqc-report-<date>/` — HTML QC reports and `images/` (read length distributions, PHRED score plots, correlation plots, etc.) used to assess sequencing quality and mapping characteristics.
  - `versions.yml` — versions of ToulligQC and dependent tools.

Notes:
- Most result folders also include a `versions.yml` file recording tool versions used during the run for reproducibility.
- Filenames include the sample name (here `nanopore_reads`) when multiple samples are processed; expect one subdirectory per sample in `minimap2/`, `salmon/`, `stringtie/`, etc.

This detailed listing should help you locate specific outputs for downstream analysis or visualization (for example, open the `execution_report_*.html` in a browser to inspect pipeline-level statistics, or load `nanopore_reads.bam` into IGV together with `combined.gtf` for visualization).

---

For more details on parameters and customization, see the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html) and the pipeline source code.
