# Sugarcane Genome Alignment using BWA-MEM

This document outlines the workflow for aligning paired-end Illumina sequencing reads from sugarcane samples against a reference genome. The process uses `bwa-mem` for alignment and `samtools` for output processing.

## Prerequisites

- **Software:** `BWA` and `Samtools` must be installed.
- **Input Files:**
    - A reference genome in FASTA format (e.g., `reference.fna`).
    - Cleaned, paired-end FASTQ files (e.g., `sample_R1.fastq.gz` and `sample_R2.fastq.gz`).
## Step 1: Reference Genome Preparation

For sugarcane genomic studies, we recommend using the high-quality, allele-defined reference genome of *Saccharum spontaneum* published by Zhang et al., 2018 in Nature Genetics.

First, download the genome from NCBI, unzip it, and then create the BWA index.

```bash
# 1.1 - Download the reference genome
wget [https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/757/435/GCF_009757435.1_ASM1313833v1/GCF_009757435.1_ASM1313833v1_genomic.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/757/435/GCF_009757435.1_ASM1313833v1/GCF_009757435.1_ASM1313833v1_genomic.fna.gz)

# 1.2 - Unzip the file
gunzip GCF_009757435.1_ASM1313833v1_genomic.fna.gz

# 1.3 - Index the reference genome (this step is slow)
bwa index GCF_009757435.1_ASM1313833v1_genomic.fna
```
## Step 2: Indexing the Reference Genome

First, an index of the reference genome must be created. This step is only required once per reference genome.

```bash
bwa index /path/to/your/reference_genome.fna
```
## Step 3: Read Alignment
This command aligns the paired-end reads to the indexed reference and pipes the output directly to samtools to generate a compressed BAM file.
```bash
bwa mem -t 8 \
-R '@RG\tID:sample_name\tSM:sample_name\tPL:ILLUMINA' \
/path/to/your/reference_genome.fna \
/path/to/your/reads_R1.fastq.gz \
/path/to/your/reads_R2.fastq.gz \
| samtools view -bS -o output_sample.bam
```


