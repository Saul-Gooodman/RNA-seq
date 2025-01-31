#!/bin/bash
#SBATCH --job-name=FastQC_job          # Job name
#SBATCH --cpus-per-task=4              # Number of CPU cores per task
#SBATCH --mem=1000                     # Memory requirement (MB)
#SBATCH --time=2:00:00                 # Time limit for the job
#SBATCH --partition=pibu_el8           # Partition

# Email notifications
#SBATCH --mail-user=yuwei.liu1@students.unibe.ch     
#SBATCH --mail-type=BEGIN,END,FAIL

# Run featureCounts using the specified Apptainer container
apptainer exec /containers/apptainer/subread_2.0.1--hed695b0_0.sif featureCounts -p -s 2 -Q 10 -t exon -g gene_id -a /data/users/yliu2/RNAseq_analysis/hisat2_mapping/Mus_musculus.GRCm39.113.gtf -o /data/users/yliu2/RNAseq_analysis/hisat2_mapping/counts_matrix.txt /data/users/yliu2/RNAseq_analysis/hisat2_mapping/*_sorted.bam
