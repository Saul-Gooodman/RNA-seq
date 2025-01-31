#!/bin/bash
#SBATCH --job-name=FastQC_job          # Job name
#SBATCH --cpus-per-task=4              # Number of CPU cores per task
#SBATCH --mem=25000                    # Memory requirement (MB)
#SBATCH --time=24:00:00                # Time limit for the job
#SBATCH --partition=pibu_el8           # Partition

# Email notifications
#SBATCH --mail-user=yuwei.liu1@students.unibe.ch     
#SBATCH --mail-type=BEGIN,END,FAIL

# Set variables
REF_GENOME="/data/users/yliu2/RNAseq_analysis/hisat2_mapping/Mus_musculus.GRCm39.dna.primary_assembly.fa"  # Path to reference genome file
GTF_FILE="/data/users/yliu2/RNAseq_analysis/hisat2_mapping/Mus_musculus.GRCm39.113.gtf"  # Path to the GTF file
INDEX_BASE="/data/users/yliu2/RNAseq_analysis/hisat2_indices/Mus_musculus_genome_index"  # Path to index files

# Mapping for each sample
for FASTQ_1 in /data/courses/rnaseq_course/toxoplasma_de/reads/*_1.fastq.gz
do
    # Get the sample name (remove the trailing "_1.fastq.gz")
    SAMPLE=$(basename $FASTQ_1 _1.fastq.gz)
    
    # Set the corresponding reverse file path
    FASTQ_2="/data/courses/rnaseq_course/toxoplasma_de/reads/${SAMPLE}_2.fastq.gz"
    
    # Set the output file paths
    OUTPUT_SAM="/data/users/yliu2/RNAseq_analysis/hisat2_mapping/${SAMPLE}_output.sam"
    OUTPUT_BAM="/data/users/yliu2/RNAseq_analysis/hisat2_mapping/${SAMPLE}_output.bam"
    SORTED_BAM="/data/users/yliu2/RNAseq_analysis/hisat2_mapping/${SAMPLE}_sorted.bam"
    THREADS=4  

    # Step 1: Run Hisat2 to map the reads to the reference genome with RF strandedness
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
      hisat2 -x $INDEX_BASE -1 $FASTQ_1 -2 $FASTQ_2 -S $OUTPUT_SAM -p $THREADS --rna-strandness RF

    # Step 2: Convert SAM to BAM format using Samtools
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
      samtools view -hbS $OUTPUT_SAM > $OUTPUT_BAM

    # Step 3: Sort the BAM file by genomic coordinates using Samtools
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
      samtools sort -m 2G -@ $THREADS -o $SORTED_BAM -T temp $OUTPUT_BAM

    # Step 4: Index the sorted BAM file using Samtools
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
      samtools index $SORTED_BAM
done