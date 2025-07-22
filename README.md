# Sugarcane Genome Alignment using BWA-MEM

This document outlines the workflow for aligning paired-end Illumina sequencing reads from sugarcane samples against a reference genome. The process uses `bwa-mem` for alignment and `samtools` for output processing.

## Prerequisites

- **Software:** `BWA` and `Samtools` must be installed.
- **Input Files:**
    - A reference genome in FASTA format (e.g., `reference.fna`).
    - Cleaned, paired-end FASTQ files (e.g., `sample_R1.fastq.gz` and `sample_R2.fastq.gz`).

## Step 1: Indexing the Reference Genome

First, an index of the reference genome must be created. This step is only required once per reference genome.

```bash
bwa index /path/to/your/reference_genome.fna

bwa mem -t 8 \
-R '@RG\tID:sample_name\tSM:sample_name\tPL:ILLUMINA' \
/path/to/your/reference_genome.fna \
/path/to/your/reads_R1.fastq.gz \
/path/to/your/reads_R2.fastq.gz \
| samtools view -bS -o output_sample.bam
