#!/bin/bash
#SBATCH --job-name=FastQC_job          # Job name
#SBATCH --cpus-per-task=1              # Number of CPU cores per task
#SBATCH --mem=1000                     # Memory requirement (MB)
#SBATCH --time=3:00:00                 # Time limit for the job
#SBATCH --partition=pibu_el8           # Partition

# Email notifications
#SBATCH --mail-user=yuwei.liu1@students.unibe.ch     
#SBATCH --mail-type=BEGIN,END,FAIL             

# Run Hisat2 to build the genome index
apptainer exec /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2-build Mus_musculus.GRCm39.dna.primary_assembly.fa Mus_musculus_genome_index