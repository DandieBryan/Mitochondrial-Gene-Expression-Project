# Mitochondrial-Gene-Expression-Project

This repository contains scripts for a pipeline designed to analyze RNA-Seq data with a focus on mitochondrial gene expression. The pipeline utilizes tools like `RNA-SeQC` and `Salmon` to calculate gene expression levels and includes steps for quality control, data compilation, and comparative analysis. The data used in this pipeline comes from lymphoblastoid cell lines provided by the [1000 Genomes Project].

## Overview

The pipeline is divided into four key components:

1. **Quality Control and Gene Expression Calculation**
   - Uses `RNA-SeQC` to perform quality control analysis on RNA-Seq data.
   - Calculates gene expression in TPM (Transcripts Per Million).

2. **Data Compilation**
   - Combines gene expression information from all samples into a single consolidated file.

3. **Comparative Analysis of Mitochondrial Gene Expression**
   - Calls an R script to compare mitochondrial gene expression calculated by `RNA-SeQC` and `Salmon`.

4. **Visualization of Gene Expression Comparisons**
   - Compares and visualizes the gene expression of mitochondrial genes using results from `RNA-SeQC` and `Salmon`.

## Prerequisites

### Software and Tools

- **RNA-SeQC**: For RNA-Seq quality control and gene expression calculation.
- **Salmon**: For quantification of transcript abundance (Transcript abundances were calculated by Salmon separately outside of this pipeline).
- **R**: For statistical analysis and visualization (requires libraries such as `ggplot2`, `dplyr`, etc.).

### Input Data

- Raw RNA-Seq data in FASTQ format from lymphoblastoid cell lines provided by the [1000 Genomes Project].
- Gene expression outputs generated by `RNA-SeQC` and `Salmon`.

## Usage

1. **Quality Control and Expression Calculation**
   - Run the first script to analyze RNA-Seq data with `RNA-SeQC` and calculate TPM values for each sample.

2. **Compile Gene Expression Data**
   - Use the second script to merge gene expression data from all samples into a single file for downstream analysis.

3. **Compare Mitochondrial Gene Expression**
   - Execute the third script to invoke the R script for comparing mitochondrial gene expression levels calculated by `RNA-SeQC` and `Salmon`.

4. **Visualize Comparative Results**
   - Run the final script to generate visualizations of the comparative analysis, highlighting differences in mitochondrial gene expression between the two tools.

## Output

The pipeline generates:

- Quality-controlled RNA-Seq data.
- TPM values for gene expression across all samples.
- Comparative metrics and statistical analysis of mitochondrial gene expression.
- Visualizations (e.g. scatter plots) comparing results from `RNA-SeQC` and `Salmon`.
