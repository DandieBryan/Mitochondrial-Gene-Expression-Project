#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --ntasks=8
#SBATCH --mem=32g
#SBATCH --tmp=32g



##########################################################################################
### MADE BY BRYAN LE, JULY 2024

# This part of the pipeline calls an R script to compare gene expression of mitochondrial 
# genes that were calculated by 'RNA-SeQC' and 'Salmon'



##########################################################################################
# Load Modules
module load R/4.4.0-openblas-rocky8

# Set paths to the input files and output plot
rnaseqc_GE_FILE="INSERT PATH TO COMPILED GE FILE GENERATED FROM 'RNA-SeQC'"
salmon_GE_FILE="INSERT PATH TO COMPILED GE FILE GENERATED FROM 'Salmon'"
output_DIRECTORY="INSERT PATH TO OUTPUT DIRECTORY"
reference_genes_FILE="INSERT PATH TO 'MT_genes.csv' FILE"

# Run the R script to generate the scatterplot
Rscript 4-MT_GE_comparison.R "$rnaseqc_GE_FILE" "$salmon_GE_FILE" "$output_DIRECTORY" "$reference_genes_FILE"