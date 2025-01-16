#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=le000209@umn.edu
#SBATCH -p moria



##########################################################################################
### Made BY BRYAN LE, JULY 2024

# The first part of this pipeline uses 'RNA-SeQC' to perform quality control analysis on
# RNA-SEQ data and calculates gene expression in TPM



##########################################################################################
# Define the directory containing the BAM files and the GTF file
rnaseq_DIRECTORY="INSERT PATH TO RNA-SEQ BAM FILES"
gtf_FILE="INSERT PATH TO GTF FILE"
rnaseqc_PATH="INSERT PATH TO 'RNA-SeQC' EXECUTABLE"
parallel_PATH="INSERT PATH TO 'Parallel' EXECUTABLE"
OUTPUT_DIR="INSERT PATH TO OUTPUT DIRECTORY"

# Change to the directory where rnaseQC output will be stored
cd $OUTPUT_DIR

# Function to run rnaseqc on a single bam file
run_rnaseqc() {
    local bam_FILE=$1
    $rnaseqc_PATH "$gtf_FILE" "$bam_FILE" --coverage .

export -f run_rnaseqc
export gtf_FILE
export rnaseqc_PATH

# Use parallel to run rnaseqc in parallel for all BAM files in the directory
find "$rnaseq_DIRECTORY" -name "*.bam" | "$parallel_PATH" run_rnaseqc

# Compile gene expression information for all samples together
sbatch 2-compile.sh