# Sugarcane Genome Alignment using BWA-MEM

Workflow for aligning paired-end Illumina sequencing reads.
The process uses `bwa-mem` for alignment and `samtools` for output processing.

## Prerequisites

- **Software:** [BWA](https://github.com/lh3/bwa.git) and [Samtools](https://github.com/samtools/samtools.git) must be installed.
- **Input Files:**
    - A reference genome in FASTA format (e.g., `reference.fna`).
    - Cleaned, paired-end FASTQ files (e.g., `sample_R1.fastq.gz` and `sample_R2.fastq.gz`).
## Step 1: Reference Genome Preparation

For sugarcane genomic studies, i recommend using the reference genome of  [Saccharum officinarum x spontaneum](https://phytozome-next.jgi.doe.gov/info/SofficinarumxspontaneumR570_v2_1) published by [Zhang et al., 2018 in Nature Genetics.](https://pmc.ncbi.nlm.nih.gov/articles/PMC11041754/)

First, download the genome from phytozome and unzip it.

```bash
# 1.1 - Download the reference genome
wget https://files.phytozome-next.jgi.doe.gov/info/SofficinarumxspontaneumR570_v2.1/v2.1/assembly/SofficinarumxspontaneumR570_v2.1.assembly.fasta.gz

# 1.2 - Unzip the file
gunzip SofficinarumxspontaneumR570_v2.1.assembly.fasta.gz

```
## Step 2: Indexing the Reference Genome

Before alignment, the reference genome FASTA file must be indexed.

```bash
bwa index /path/to/your/reference_genome.fna
```
## Step 3: Read Alignment
This command aligns the paired-end reads to the indexed reference and the output directly to samtools to generate a compressed BAM file.
```bash
bwa mem -t 8 \
-R '@RG\tID:sample_name\tSM:sample_name\tPL:ILLUMINA' \
/path/to/your/reference_genome.fna \
/path/to/your/reads_R1.fastq.gz \
/path/to/your/reads_R2.fastq.gz \
| samtools view -bS -o output_sample.bam
```


