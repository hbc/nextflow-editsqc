# Nextflow Nanopore Pipeline

## Introduction

This pipeline provides a reproducible workflow for processing Oxford Nanopore sequencing data using Nextflow and nf-core modules. It automates key steps from raw read filtering to transcript quantification and report generation, leveraging containers for portability across systems.

It is aimed at small genomes with optional plasmids: the genome and plasmid sequences and annotations are combined into a single reference, reads are aligned to both the genome and the transcriptome, and transcripts are assembled and quantified against the reference annotation.

## Requirements

- [Nextflow](https://www.nextflow.io/) `>=24.04`
- Java 17 or later
- [Docker](https://www.docker.com/) or Singularity/Apptainer

All tools run inside containers, so nothing else needs to be installed locally. Every example below uses `-profile docker`; substitute `singularity` if that is what your system provides.

## Pipeline Summary

The pipeline performs the following steps:

1. **Reference preparation**
   Converts the genome and (optional) plasmid GFF annotations to GTF, concatenates the genome and plasmid FASTA and GTF files into a single combined reference, then sorts, bgzips and indexes it (`gffread`, `samtools faidx`, `tabix`).

2. **Transcript extraction**
   Derives a transcript-level GTF and a transcriptome FASTA from the combined reference (`gffread`).

3. **Read filtering**
   Filters raw reads by length (`fastplong`, default `--length_limit 10000`).

4. **Quality control of raw reads**
   Runs ToulligQC on the filtered reads.

5. **Alignment**
   Aligns reads to both the combined genome and the transcriptome with `minimap2`, then indexes the genome alignments and collects `idxstats`.

6. **Region reports** *(optional, per sample)*
   For samples with a `region` value, subsets the genome BAM to that region, builds a BED of the features it contains, and generates a self-contained IGV HTML report.

7. **Quality control of mapped reads**
   Keeps only primary mapped reads with MAPQ >= 50 and runs ToulligQC on them.

8. **Transcriptome assembly and quantification**
   Assembles transcripts with StringTie, merges them with the reference annotation, restores reference gene IDs, and re-quantifies against the merged annotation to produce final abundances.

## Supplying Inputs

Inputs are supplied through a samplesheet CSV, passed with `--input`:

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results
```

### Samplesheet format

One row per sample, with this header:

```csv
sample,fastq_files,genome_fasta,gene_annotation,plasmid_fasta,plasmid_annotation,region
sample1,reads/sample1.fastq.gz,ref/genome.fasta,ref/genome.gff,ref/plasmid.fasta,ref/plasmid.gff,plasmid
sample2,reads/sample2.fastq.gz,ref/genome.fasta,ref/genome.gff,,,
```

| Column | Required | Description |
| --- | --- | --- |
| `sample` | yes | Sample identifier. Used to name every output file for that sample, so keep it short and free of spaces. |
| `fastq_files` | yes | Path to the FASTQ file for this sample. |
| `genome_fasta` | yes | Reference genome FASTA. |
| `gene_annotation` | yes | Genome annotation in GFF format. |
| `plasmid_fasta` | no | Plasmid FASTA. Leave blank to run genome-only. |
| `plasmid_annotation` | no | Plasmid annotation in GFF format. Leave blank to run genome-only. |
| `region` | no | Sequence name to report on, e.g. `plasmid`. Leave blank to skip the region and IGV report steps for that sample. |

Reference files may differ between rows, so samples with different genomes or plasmids can be processed in the same run.

`region` is evaluated per sample: in the example above `sample1` produces a region BAM and an IGV report for `plasmid`, while `sample2` skips those steps entirely. The value must match a sequence name present in the combined reference for that sample — a genome-only sample cannot report on `plasmid`.

### Other parameters

- `--outdir` : Output directory. Default: `./results`
- `--valid_features` : Space-separated list of feature types in the GFF to keep during conversion. Default: `"gene CDS transcription_unit tRNA rRNA"`
- `--extension` : Number of bases to extend features by during GFF-to-GTF conversion. Default: `0`
- `--single_end` : Whether reads are single-end. Default: `true`

Parameters can also be collected in a config file:

```groovy
params {
    input          = 'samplesheet.csv'
    outdir         = './results'
    valid_features = 'gene CDS tRNA rRNA'
    extension      = '0'
}
```

```bash
nextflow run main.nf -profile docker -c mysample.config
```

### Specific details about FASTA and GFF

Below are specific formatting rules and examples to ensure your FASTA and GFF inputs are compatible with the pipeline.

- FASTA (genome or transcriptome)
  - Headers must start with a single '>' character followed by a unique sequence identifier (sequence IDs are used to match annotation seqids). Avoid spaces in the primary ID; use underscores if needed.
    - Good: `>chromosome` or `>plasmid_1` or `>TAC.gene01.1`
    - Avoid: `>chromosome 1 description` (the space makes the full header the ID)
  - The sequence identifier in the FASTA header must exactly match the seqid fields used in your GFF (case-sensitive).
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
  - The `region` column in the samplesheet must also match a sequence name in the combined reference, or the region report will come back empty.
  - Ensure no duplicate sequence IDs in FASTA; tools (StringTie) require unique IDs.
  - GFF3 multi-line attributes or unusual quoting can sometimes break parsers; provide standard, single-line attribute fields.
  - Provide a GFF3 — if you have GFF3 the pipeline will convert to GTF for StringTie where needed.

## Running the pipeline

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results
```

### Test profiles

Two profiles run the pipeline on the small dataset in `data/` and are the quickest way to check an installation:

```bash
# genome plus plasmid, with a region report
nextflow run main.nf -profile test,docker

# genome only, blank region: the region and IGV steps are skipped
nextflow run main.nf -profile test_genome,docker
```

Both cap CPU, memory and runtime in `conf/test.config` and `conf/test_genome.config` so they fit on a laptop or a CI runner.

### Other profiles

- `docker`, and the container settings for Singularity/Apptainer, are defined in `nextflow.config`.
- `minimap` switches genome alignment to spliced mode (`-ax splice -uf -k14`) instead of the default `-ax map-ont`.

## Pipeline Outputs

All results are written to `--outdir` (default `results/`). Output files are named after the `sample` column of the samplesheet, so a sample called `sample1` produces `sample1.bam`, `sample1.gene.abundance.txt`, and so on. Reference-derived files are named after the input FASTA, for example `genome_plasmid_combined.fasta`.

```text
results/
├── combine_fasta_annotation/   combined reference FASTA, GTF and indexes
├── filtered_reads/             length-filtered FASTQ
├── gffread/                    transcript GTF and transcriptome FASTA
├── igvreports/                 per-sample IGV HTML reports (region samples only)
├── minimap2/
│   ├── genome/                 genome alignments, idxstats, QC, region subsets
│   └── tx/                     transcriptome alignments
├── stringtie/                  initial assembly and merged annotation
├── stringtie_new/              final quantification
├── toulligqc/                  QC of filtered reads
└── pipeline_info/              execution reports, timeline, trace and DAG
```

### Most useful files

- `stringtie_new/`
  - `<sample>.gene.abundance.txt` — gene-level abundance table with Coverage, FPKM and TPM columns. This is the main quantification result.
  - `<sample>.transcripts.gtf` — transcripts quantified against the merged annotation.
  - `<sample>.coverage.gtf` — per-locus coverage.
  - `<sample>.ballgown/` — Ballgown-compatible tables for downstream expression analysis.

- `minimap2/genome/`
  - `<sample>.bam` and `<sample>.bam.bai` — reads aligned to the combined reference.
  - `<sample>.idxstats` — mapped and unmapped read counts per sequence, the quickest way to see how reads split between genome and plasmid.
  - `qc/` — ToulligQC report for primary mapped reads (MAPQ >= 50).
  - `region/<sample>.<region>.bam` — alignments subset to the requested region, for samples that specify one.

- `combine_fasta_annotation/`
  - `<genome>_<plasmid>_combined.fasta` — combined reference used for alignment, plus its `.fai` index.
  - `<genome>_<plasmid>_combined.gtf` — combined annotation, plus a sorted, bgzipped and tabix-indexed copy for use as an IGV track.
  - `regions.bed` — features contained in the requested region.

- `igvreports/`
  - `<sample>_report.html` — self-contained IGV report for the requested region, viewable in a browser with no additional files.

### Complementary output files

- `filtered_reads/`
  - `<reads>_filtered.fastq.gz` — reads after length filtering by `fastplong`.

- `gffread/`
  - `<reference>_fixed.gtf` — transcript-level GTF derived from the combined annotation, used by StringTie.
  - `<reference>_tx.fasta` — transcriptome FASTA used for transcriptome alignment.

- `minimap2/tx/`
  - `<sample>.bam` and `<sample>.bam.bai` — reads aligned to the transcriptome.

- `stringtie/`
  - `<sample>.transcripts.gtf`, `<sample>.coverage.gtf`, `<sample>.gene.abundance.txt` — results of the initial assembly, before merging.
  - `stringtie.merged.gtf` — assembled transcripts merged with the reference annotation.
  - `stringtie.merged.refgene.gtf` — the merged annotation with reference gene IDs restored. This is what the final quantification runs against.

- `toulligqc/`
  - `<sample>Toulligqc-report-<date>/` — HTML QC report and images (read length distributions, PHRED score plots) for the filtered reads.

- `pipeline_info/`
  - `execution_report_<timestamp>.html` — runtime, resource usage and process summaries.
  - `execution_timeline_<timestamp>.html` — timeline of workflow execution.
  - `execution_trace_<timestamp>.txt` — raw task trace.
  - `pipeline_dag_<timestamp>.html` — graphical representation of the pipeline DAG.

Most result folders also include a `versions.yml` recording the tool versions used, for reproducibility.

To inspect a run, open `pipeline_info/execution_report_*.html` in a browser, or load `minimap2/genome/<sample>.bam` into IGV together with `combine_fasta_annotation/<reference>_combined.gtf`.

## Credits

This pipeline was written by:

- [Lorena Pantano](https://www.linkedin.com/in/lpantano/)
- [Alex Bartlett](https://www.linkedin.com/in/alexandra-bartlett-926b32109/)

In collaboration with:

- Alkmini Diamantidi (<alkmini.diamantidi@braskem.com>)
- Hugo Rody Vianna Silva (<hugo.vianna@braskem.com>)
- Susan McAvoy (<susan.mcavoy@braskem.com>)
- Jonathan Turner (<jonjoet@gmail.com>)
- Wing-On Ng (<wingon.ng@braskem.com>)

## Citations

Tool references are collected in [CITATIONS.md](CITATIONS.md). Release history is in [CHANGELOG.md](CHANGELOG.md), and the pipeline is released under the [MIT license](LICENSE).

---

For more details on parameters and customization, see the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html) and the pipeline source code.
