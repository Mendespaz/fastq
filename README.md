# fastq

## Alinhamento 
bwa mem -t 8 -R '@RG\tID:pai_alto\tSM:pai_alto\tPL:ILLUMINA' Cana/genoma_R570/GCA_038087645.1_Saccharum_officinarum_X_spontaneum_var_R570_v2.1_genomic.fna Cana/pais/cat/cleaned_reds/pai_alto_R1_cleaned_paired.fastq.gz Cana/pais/cat/cleaned_reds/pai_alto_R2_cleaned_paired.fastq.gz | samtools view -bS -o pai_alto.bam
