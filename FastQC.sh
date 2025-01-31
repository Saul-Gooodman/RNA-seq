#!/bin/bash
#SBATCH --job-name=FastQC_job          # Job name
#SBATCH --cpus-per-task=1              # Number of CPU cores per task
#SBATCH --mem=1000                     # Memory requirement (MB)
#SBATCH --time=3:00:00                 # Time limit for the job
#SBATCH --partition=pibu_el8           # Partition

# Email notifications
#SBATCH --mail-user=yuwei.liu1@students.unibe.ch     
#SBATCH --mail-type=BEGIN,END,FAIL             

# Run FastQC container for quality control
apptainer exec --bind /data/ /containers/apptainer/fastqc-0.12.1.sif fastqc -o /data/users/yliu2/RNAseq_analysis/fastqc_results /data/courses/rnaseq_course/toxoplasma_de/reads/*.fastq.gz