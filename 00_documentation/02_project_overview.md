# SMART-Seq2 Data Analysis Pipeline - Project Overview

## Project Information

| Field | Value |
|-------|-------|
| **Project Name** | End-to-End Single Cell RNA Sequencing Data Analysis Pipeline (SMART-Seq2) |
| **Author** | Sourav Roy (s2599932) |
| **Institution** | University of Edinburgh, MSc Data Science for Biology |
| **Repository Status** | Git repository, main branch |

---

## Project Purpose and Scope

This MSc dissertation project implements a complete single-cell RNA sequencing (scRNA-seq) analysis pipeline using the SMART-Seq2 protocol. The project addresses three main objectives:

1. **Replication Phase:** Re-evaluate and replicate findings from previous mouse studies at the Tomlinson Laboratory regarding stably expressed genes (SEGs)

2. **Upstream Processing:** Analyze raw FASTQ files from human single-cell RNA sequencing datasets using SMART-Seq2 protocol

3. **Downstream Analysis:** Implement a refined algorithmic approach to identify stably expressed genes (SEGs) that differ between cell types/developmental stages, particularly examining transcriptional changes between human embryonic stem cells (hESCs) and early human development cells

---

## Directory Structure

```
SMART-Seq2_Data_Analysis_Pipeline/
├── .git/                              # Git repository files
├── .gitignore                         # Git ignore rules
├── README.md                          # Project README
│
├── 00_documentation/                  # Project documentation
│   ├── 00_planned_changes.md         # Plans for project reorganization
│   ├── 01_background.md              # Background and methodology
│   ├── 02_method                     # Method details (stub)
│   └── 03_project_overview.md        # This file
│
├── 01_enviroment-setup/               # Environment configuration
│   ├── conda/
│   │   ├── enviroment_upstream.yaml  # Upstream analysis conda environment
│   │   └── enviroment_downstream.yaml # Downstream analysis conda environment
│   └── dissertation-sessioninfo-files/
│       └── finalized-human-ranked-SEG-enrichment_session_info.txt
│
├── 02_data/                           # Data directory (sparse by design)
│   └── .gitignore                    # Ignores large data files (.rds, .RData)
│
└── 03_analysis-pipelines/             # Core analysis code
    ├── 01_upstream-analysis/          # Raw reads to counts matrix
    │   ├── pre-processing-pipeline.Rmd
    │   └── metadata-creation/
    │       ├── 01-targets-file-creation.Rmd
    │       ├── 02-metadata-from-sdrf.Rmd
    │       └── GSMtoSRR/
    │           └── GSM-SRR.sh
    │
    ├── 02_downstream-analysis/        # Analysis of processed data
    │   ├── revised-human-seg-algorithm.R
    │   └── pouf1.R
    │
    └── miscallenous-analysis/         # Exploratory & supplementary
        ├── investigate_pca_discrepancy.R
        ├── reproduce_findings.R
        └── figures/
            ├── QC/
            └── other_plots/
```

---

## Technologies and Tools

### Bioinformatics Tools (Upstream Processing)

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | 0.11.9 | Read quality assessment |
| MultiQC | 1.22.3 | Quality report aggregation |
| Trimmomatic | 0.39 | Read trimming and quality filtering |
| STAR | 2.7.10a | RNA-seq read alignment |
| featureCounts | 2.0.3 | Transcript quantification |
| samtools | - | BAM file manipulation and QC |
| bedtools | - | Genomic analysis utilities |

### R/Bioconductor Ecosystem (Downstream Analysis)

**Core Packages:**
- **Seurat** (v5.1.0): Single-cell analysis framework, clustering, dimensionality reduction
- **SingleCellExperiment/SummarizedExperiment**: Data structure management
- **scran** (v1.32.0): Single-cell RNA-seq analysis (normalization, QC)
- **scater** (v1.32.0): Single-cell analysis, QC visualization
- **scMerge** (v1.20.0): Source for stably expressed genes (SEG) reference lists

**Specialized Analysis:**
- **fgsea** (v1.30.0): Fast gene set enrichment analysis
- **msigdb** (v1.12.0): Molecular signatures database
- **biomaRt** (v2.60.1): Gene annotation and mapping

**Visualization:**
- **ggplot2** (v3.5.1): Publication-quality graphics
- **pheatmap** (v1.0.12): Heatmaps
- **cowplot** (v1.1.3): Multi-panel plots

**Data Manipulation:**
- **tidyverse** ecosystem: dplyr, tidyr, readr, stringr, purrr
- **data.table** (v1.15.4): Fast data frame operations

### Environment Management
- **Conda** (v4.4.1): Package management
- **R** (v4.4.1): Statistical computing

---

## Analysis Workflows

### 1. Upstream Analysis (FASTQ to Counts Matrix)

**File:** `03_analysis-pipelines/01_upstream-analysis/pre-processing-pipeline.Rmd`

**Pipeline Steps:**
1. **Quality Assessment:** FastQC + MultiQC
2. **Trimming:** Trimmomatic (removes low-quality bases, crop to 50bp)
3. **Reference Preparation:** Download human genome (GRCh38) and annotations from Ensembl
4. **Alignment:** STAR (maps reads to reference genome)
5. **Quality Check:** samtools flagstat (alignment statistics)
6. **Quantification:** featureCounts (counts reads mapping to genes)
7. **Post-processing:** Clean column names and format output

**Output:** Gene expression count matrix (genes × samples)

### 2. Metadata Creation

**Files:**
- `01-targets-file-creation.Rmd`: Creates sample metadata from GEO
- `02-metadata-from-sdrf.Rmd`: Alternative metadata extraction from SDRF
- `GSM-SRR.sh`: Maps GEO sample IDs (GSM) to SRA run IDs (SRR)

### 3. Downstream Analysis (Gene Expression Analysis)

**Main File:** `03_analysis-pipelines/02_downstream-analysis/revised-human-seg-algorithm.R`

**Workflow:**
1. **Data Loading:** Load count matrices and create SingleCellExperiment objects
2. **Quality Control:** Filter cells based on library size, gene detection, mitochondrial content
3. **Normalization:** Deconvolution normalization using scran
4. **Dimensionality Reduction:** PCA via Seurat
5. **Clustering:** Leiden algorithm clustering
6. **Gene Expression Analysis:** Identify differentially stable genes between conditions
7. **Enrichment Analysis:** FGSEA pathway analysis
8. **Visualization:** PCA plots, violin plots, heatmaps

### 4. Supplementary Analysis

| Script | Purpose |
|--------|---------|
| `investigate_pca_discrepancy.R` | Compare PCA implementations between packages |
| `reproduce_findings.R` | Replicate original mouse study findings |
| `pouf1.R` | Gene-specific analysis for POU5F1 (Oct4) |

---

## Key Algorithm: SEG-based Cell Type Differentiation

The novel algorithmic approach in this project:

1. Identifies genes that are **stably expressed** (low variance) in one cell type but **highly variable** in another
2. Uses ranked gene lists to perform functional enrichment analysis (FGSEA)
3. Reveals transcriptomic differences at the pathway level
4. Combines concepts from Yingxin et al. 2019 with modern enrichment analysis

### Quality Control Standards
- 3 MAD (median absolute deviations) as outlier threshold
- Filters: total counts, detected genes, mitochondrial percentage
- ERCC spike-in controls for mouse data

---

## Datasets Analyzed

| Dataset | Type | Description |
|---------|------|-------------|
| E-MTAB-2600 | Mouse | Embryonic stem cells (replication/validation) |
| E-MTAB-3929 | Human | Early development stages (main analysis) |
| GSE75748 | Human | hESC SMART-Seq2 data |
| GSE36552 | Human | Reference dataset |
| GSE71318 | Human | Reference dataset |

---

## Conda Environments

### Upstream Environment (`enviroment_upstream.yaml`)
- **Purpose:** Raw read processing and quality control
- **Key Tools:** FastQC, MultiQC, Trimmomatic, STAR, featureCounts, samtools
- **R Packages:** tidyverse, biomaRt, GEOquery

### Downstream Environment (`enviroment_downstream.yaml`)
- **Purpose:** Single-cell data analysis and visualization
- **Key Tools:** Seurat, scran, scater, fgsea
- **Specialized:** scMerge (SEG references), ComplexHeatmap

---

## Output Figures

### Quality Control Plots (`figures/QC/`)
- `total_counts.png` - Library size distribution
- `detected_genes.png` - Gene detection per cell
- `mitochondrial_proportion.png` - Mitochondrial content
- `ERCC_proportion.png` - Spike-in proportion

### Analysis Plots (`figures/other_plots/`)
- `all_samples_pca.png` - PCA of all samples
- `cluster_wise_pca.png` - PCA by cluster
- `nono_violin.png` - NONO gene expression
- `seg_cv.png` - SEG coefficient of variation

---

## File Summary

| File Type | Count | Purpose |
|-----------|-------|---------|
| R Markdown (.Rmd) | 3 | Documented workflows |
| R Scripts (.R) | 4 | Analysis code |
| Shell Scripts (.sh) | 1 | ID mapping utility |
| YAML | 2 | Conda environments |
| Markdown | 4 | Documentation |
| PNG Images | 8 | Visualizations |

---

## Reproducibility Notes

1. **Environment Recreation:**
   ```bash
   # Upstream analysis
   conda env create -f 01_enviroment-setup/conda/enviroment_upstream.yaml

   # Downstream analysis
   conda env create -f 01_enviroment-setup/conda/enviroment_downstream.yaml
   ```

2. **Session Info:** R package versions recorded in `dissertation-sessioninfo-files/`

3. **Data:** Users must populate `02_data/` with their own datasets (large files excluded from git)

---

## Project Status

**Current State:** Analysis complete for both mouse (replication) and human (original) datasets

**Planned Improvements:** (see `00_planned_changes.md`)
- Rewrite upstream pipeline in Nextflow for HPC compatibility
- Convert downstream analysis to modular R Markdown format
- Comprehensive documentation updates

---

*Document generated: 2026-02-04*