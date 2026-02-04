# SMART-Seq2 Single-Cell RNA Sequencing Analysis Pipeline

A comprehensive end-to-end pipeline for analyzing single-cell RNA sequencing data using the SMART-Seq2 protocol. This project implements both upstream processing (raw FASTQ to count matrix) and downstream analysis (quality control, normalization, clustering, and gene set enrichment analysis).

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Project Structure](#project-structure)
- [Usage](#usage)
  - [Upstream Analysis](#upstream-analysis)
  - [Downstream Analysis](#downstream-analysis)
- [Datasets](#datasets)
- [Output](#output)
- [Documentation](#documentation)
- [Author](#author)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Overview

This pipeline was developed as part of an MSc dissertation project at the University of Edinburgh to:

1. **Replicate** findings from previous mouse studies on stably expressed genes (SEGs)
2. **Process** raw FASTQ files from human SMART-Seq2 experiments
3. **Identify** stably expressed genes that differentiate between cell types and developmental stages

The novel algorithmic approach identifies genes that are stably expressed in one cell type but highly variable in another, using ranked gene lists and functional enrichment analysis to reveal transcriptomic differences at the pathway level.

## Features

- **Complete Upstream Pipeline**
  - Quality assessment with FastQC and MultiQC
  - Read trimming with Trimmomatic
  - Alignment to reference genome using STAR
  - Transcript quantification with featureCounts

- **Comprehensive Downstream Analysis**
  - Quality control with outlier detection (MAD-based filtering)
  - Normalization using scran deconvolution
  - Dimensionality reduction (PCA) and clustering (Leiden algorithm)
  - Stably expressed gene (SEG) identification
  - Gene set enrichment analysis (FGSEA)

- **Reproducible Environment**
  - Conda environment specifications for both analysis phases
  - Session info tracking for R package versions
  - Modular, documented workflows

## Prerequisites

- Linux operating system (tested on Ubuntu)
- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/)
- At least 32GB RAM (for STAR alignment)
- Sufficient storage for reference genomes and FASTQ files

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
   # Create subdirectories for your data
   mkdir -p 02_data/raw_fastq
   mkdir -p 02_data/processed
   mkdir -p 02_data/metadata
   ```

## Project Structure

```
SMART-Seq2_Data_Analysis_Pipeline/
│
├── 00_documentation/              # Project documentation
│   ├── 00_planned_changes.md     # Future improvements
│   ├── 01_background.md          # Scientific background
│   ├── 02_method                 # Methodology details
│   └── 03_project_overview.md    # Comprehensive project overview
│
├── 01_enviroment-setup/           # Environment configuration
│   ├── conda/
│   │   ├── enviroment_upstream.yaml    # Upstream tools
│   │   └── enviroment_downstream.yaml  # R/Bioconductor packages
│   └── dissertation-sessioninfo-files/ # R session info
│
├── 02_data/                       # Data directory (user-populated)
│   └── .gitignore                # Excludes large data files
│
├── 03_analysis-pipelines/         # Analysis code
│   ├── 01_upstream-analysis/     # FASTQ processing
│   │   ├── pre-processing-pipeline.Rmd
│   │   └── metadata-creation/    # Sample metadata scripts
│   │
│   ├── 02_downstream-analysis/   # Expression analysis
│   │   ├── revised-human-seg-algorithm.R
│   │   └── pouf1.R
│   │
│   └── miscallenous-analysis/    # Supplementary analyses
│       ├── investigate_pca_discrepancy.R
│       ├── reproduce_findings.R
│       └── figures/              # Output plots
│
└── README.md
```

## Usage

### Upstream Analysis

The upstream pipeline processes raw FASTQ files into a gene expression count matrix.

1. **Activate the environment**
   ```bash
   conda activate upstream_analysis
   ```

2. **Prepare metadata**

   Create a targets file with sample information:
   ```bash
   cd 03_analysis-pipelines/01_upstream-analysis/metadata-creation
   # Run the appropriate Rmd file for your data source
   ```

3. **Run the preprocessing pipeline**

   Open and execute `pre-processing-pipeline.Rmd` in RStudio or run via command line:
   ```bash
   cd 03_analysis-pipelines/01_upstream-analysis
   Rscript -e "rmarkdown::render('pre-processing-pipeline.Rmd')"
   ```

   The pipeline will:
   - Run FastQC on raw reads
   - Trim reads with Trimmomatic
   - Align to GRCh38 reference using STAR
   - Quantify transcripts with featureCounts

### Downstream Analysis

The downstream analysis performs quality control, normalization, and identifies differentially stable genes.

1. **Activate the environment**
   ```bash
   conda activate downstream_analysis
   ```

2. **Run the SEG analysis**
   ```bash
   cd 03_analysis-pipelines/02_downstream-analysis
   Rscript revised-human-seg-algorithm.R
   ```

   The analysis will:
   - Load count matrices and metadata
   - Filter low-quality cells (library size, gene detection, mitochondrial content)
   - Normalize expression data
   - Perform PCA and clustering
   - Identify stably expressed genes per cell type
   - Run gene set enrichment analysis

### Quality Control Parameters

Default QC thresholds (adjustable in the scripts):

| Metric | Threshold | Method |
|--------|-----------|--------|
| Low library size | < 3 MAD from median | Outlier detection |
| Low gene count | < 3 MAD from median | Outlier detection |
| High mitochondrial % | > 3 MAD from median | Outlier detection |

## Datasets

This pipeline has been tested with the following datasets:

| Dataset | Species | Description | Source |
|---------|---------|-------------|--------|
| E-MTAB-2600 | Mouse | Embryonic stem cells | ArrayExpress |
| E-MTAB-3929 | Human | Early development | ArrayExpress |
| GSE75748 | Human | hESC (SMART-Seq2) | GEO |
| GSE36552 | Human | Reference dataset | GEO |
| GSE71318 | Human | Reference dataset | GEO |

## Output

### Quality Control Figures
- `figures/QC/total_counts.png` - Library size distribution
- `figures/QC/detected_genes.png` - Genes detected per cell
- `figures/QC/mitochondrial_proportion.png` - Mitochondrial content
- `figures/QC/ERCC_proportion.png` - Spike-in proportion

### Analysis Figures
- `figures/other_plots/all_samples_pca.png` - PCA visualization
- `figures/other_plots/cluster_wise_pca.png` - Cluster-specific PCA
- `figures/other_plots/seg_cv.png` - SEG coefficient of variation

### Data Outputs
- Count matrix (genes x samples)
- Normalized expression matrix
- Ranked gene lists for enrichment analysis
- Enrichment analysis results (FGSEA output)

## Documentation

Detailed documentation is available in the `00_documentation/` directory:

- **Background**: Scientific context and motivation (`01_background.md`)
- **Project Overview**: Complete technical documentation (`03_project_overview.md`)
- **Planned Changes**: Future development roadmap (`00_planned_changes.md`)

## Author

**Sourav Roy**
MSc Data Science for Biology
University of Edinburgh
Student ID: s2599932

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Tomlinson Laboratory, University of Edinburgh
- scMerge package developers for stably expressed gene references (Yingxin et al., 2019)
- Bioconductor community for single-cell analysis tools

---

## Citation

If you use this pipeline in your research, please cite:

```
Roy, S. (2024). End-to-End Single Cell RNA Sequencing Data Analysis Pipeline.
MSc Dissertation, University of Edinburgh.
```

## References

- Yingxin, L., et al. (2019). scMerge leverages factor analysis, stable expression, and pseudoreplication to merge multiple single-cell RNA-seq datasets. PNAS.
- Lun, A.T., McCarthy, D.J., & Marioni, J.C. (2016). A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor. F1000Research.
- Stuart, T., et al. (2019). Comprehensive Integration of Single-Cell Data. Cell.
