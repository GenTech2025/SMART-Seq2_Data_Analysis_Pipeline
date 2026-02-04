# SMART-Seq2 Single-Cell RNA Sequencing Analysis Pipeline

**TL;DR**
An end-to-end SMART-Seq2 single-cell RNA-seq analysis pipeline with a specific focus on identifying *differentially stable genes* across cell types and developmental stages. Instead of relying solely on differential expression, this approach highlights genes that are stable in one biological context but variable in another, enabling pathway-level biological interpretation.

A comprehensive pipeline for analyzing single-cell RNA sequencing data using the SMART-Seq2 protocol. This project implements both upstream processing (raw FASTQ to count matrix) and downstream analysis (quality control, normalization, clustering, and gene set enrichment analysis).

## Table of Contents

* [Overview](#overview)
* [Key Findings](#key-findings)
* [Features](#features)
* [Prerequisites](#prerequisites)
* [Installation](#installation)
* [Project Structure](#project-structure)
* [Usage](#usage)
* [Datasets](#datasets)
* [Output](#output)
* [Documentation](#documentation)
* [Author](#author)
* [References](#references)

## Overview

This pipeline was developed as part of an MSc dissertation project at the University of Edinburgh, investigating transcriptional differences between human embryonic stem cells (hESCs) and early developmental cells using a novel stably expressed gene (SEG)–based approach.

### Project Objectives

1. **Replicate** findings from previous mouse studies on stably expressed genes (SEGs)
2. **Process** raw FASTQ files from human SMART-Seq2 experiments
3. **Identify** genes that show differential expression stability between cell types
4. **Reveal** biological pathways through gene set enrichment analysis

### The Algorithm

The approach identifies genes that are stably expressed in one cell type but highly variable in another. Genes are ranked based on differences in their stability index, and Fast Gene Set Enrichment Analysis (FGSEA) is applied to reveal transcriptomic differences at the pathway level rather than relying solely on individual gene expression changes.

## Key Findings

* Successfully validated the SEG-based approach on mouse embryonic stem cell datasets
* Pluripotency master regulators (OCT4, NANOG, SOX2) maintain stable expression in undifferentiated hESCs
* These genes become destabilized during differentiation to embryonic day 7
* Enrichment analysis revealed coordinated changes in functional gene modules during cell state transitions

See [03_results.md](00_documentation/03_results.md) for detailed results.

## Features

### Upstream Pipeline

* Quality assessment with FastQC and MultiQC
* Read trimming with Trimmomatic
* Alignment to the GRCh38 reference genome using STAR
* Transcript quantification with featureCounts

### Downstream Analysis

* Quality control with MAD-based outlier detection
* Normalization using scran deconvolution
* Stably expressed gene (SEG) identification using scMerge
* Dimensionality reduction (PCA) and visualization
* Gene set enrichment analysis (FGSEA) with MSigDB

### Reproducibility

* Conda environment specifications for both analysis phases
* R session information tracking
* Modular, documented workflows

## Prerequisites

* Linux operating system (tested on Ubuntu)
* [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/)
* At least 32 GB RAM (required for STAR alignment)
* Sufficient storage for reference genomes (~30 GB) and FASTQ files

## Installation

1. **Clone the repository**

   ```bash
   git clone https://github.com/yourusername/SMART-Seq2_Data_Analysis_Pipeline.git
   cd SMART-Seq2_Data_Analysis_Pipeline
   ```

2. **Create the upstream analysis environment**

   ```bash
   conda env create -f 01_enviroment-setup/conda/enviroment_upstream.yaml
   conda activate upstream_analysis
   ```

3. **Create the downstream analysis environment**

   ```bash
   conda env create -f 01_enviroment-setup/conda/enviroment_downstream.yaml
   conda activate downstream_analysis
   ```

4. **Prepare the data directory**

   ```bash
   mkdir -p 02_data/raw_fastq
   mkdir -p 02_data/processed
   mkdir -p 02_data/metadata
   ```

## Project Structure

```
SMART-Seq2_Data_Analysis_Pipeline/
│
├── 00_documentation/
│   ├── 01_background.md
│   ├── 02_project_overview.md
│   └── 03_results.md
│
├── 01_enviroment-setup/
│   ├── conda/
│   │   ├── enviroment_upstream.yaml
│   │   └── enviroment_downstream.yaml
│   └── dissertation-sessioninfo-files/
│
├── 02_data/
│
├── 03_analysis-pipelines/
│   ├── 01_upstream-analysis/
│   │   ├── pre-processing-pipeline.Rmd
│   │   └── metadata-creation/
│   │
│   ├── 02_downstream-analysis/
│   │   ├── revised-human-seg-algorithm.R
│   │   └── pouf1.R
│   │
│   └── miscellaneous-analysis/
│
└── README.md
```

## Usage

### Upstream Analysis

The upstream pipeline converts raw FASTQ files into a gene expression count matrix.

1. **Activate the environment**

   ```bash
   conda activate upstream_analysis
   ```

2. **Configure paths** in `pre-processing-pipeline.Rmd`

3. **Run pipeline steps**

   * FastQC + MultiQC
   * Trimmomatic
   * Reference genome download
   * STAR genome indexing and alignment
   * featureCounts quantification
   * R-based post-processing

### Downstream Analysis

1. **Activate the environment**

   ```bash
   conda activate downstream_analysis
   ```

2. **Run the SEG analysis**

   ```bash
   cd 03_analysis-pipelines/02_downstream-analysis
   Rscript revised-human-seg-algorithm.R
   ```

   This performs:

   * Quality control filtering (3 MAD threshold)
   * scran-based normalization
   * SEG index calculation
   * Differential stability analysis
   * FGSEA pathway enrichment

### Quality Control Thresholds

| Metric          | Threshold | Type        |
| --------------- | --------- | ----------- |
| Library size    | 3 MAD     | Lower bound |
| Detected genes  | 3 MAD     | Lower bound |
| Mitochondrial % | 3 MAD     | Upper bound |

## Datasets

| Dataset     | Species | Description            | Source       |
| ----------- | ------- | ---------------------- | ------------ |
| E-MTAB-2600 | Mouse   | Embryonic stem cells   | ArrayExpress |
| E-MTAB-3929 | Human   | Early development (E7) | ArrayExpress |
| GSE75748    | Human   | Undifferentiated hESCs | GEO          |

## Output

* QC and PCA visualizations
* SEG index tables and differential stability scores
* FGSEA enrichment results
* Heatmaps of stable and unstable genes

## Documentation

Detailed documentation is available in `00_documentation/`.

## Author

**Sourav Roy**
MSc Data Science for Biology
University of Edinburgh

## License

MIT License

## Acknowledgments

- [Tomlinson Laboratory](https://regenerative-medicine.ed.ac.uk/research/simon-tomlinson), University of Edinburgh
- scMerge package developers ([Yingxin et al., 2019](https://doi.org/10.1093/gigascience/giz106))
- Bioconductor community

## References

1. Lin, Y., et al. (2019). Evaluating stably expressed genes in single cells. *GigaScience*, 8(9), giz106. [DOI](https://doi.org/10.1093/gigascience/giz106)

2. Lin, Y., et al. (2019). scMerge leverages factor analysis, stable expression, and pseudoreplication to merge multiple single-cell RNA-seq datasets. *PNAS*, 116(20), 9775-9784.

3. Lun, A.T., McCarthy, D.J., & Marioni, J.C. (2016). A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor. *F1000Research*, 5, 2122.

4. Stuart, T., et al. (2019). Comprehensive Integration of Single-Cell Data. *Cell*, 177(7), 1888-1902.

---

**Citation**

If you use this pipeline in your research, please cite:

```
Roy, S. (2024). End-to-End Single Cell RNA Sequencing Data Analysis Pipeline.
MSc Dissertation, University of Edinburgh.
```
