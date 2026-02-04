# Project Results

## Executive Summary

This project successfully implemented a novel algorithm to identify transcriptional differences between cell types by analyzing stably expressed genes (SEGs). The analysis compared human embryonic stem cells (hESCs) with early developmental cells, revealing key biological pathways that differentiate pluripotent and differentiating states.

---

## Objectives and Outcomes

| Objective | Status | Key Finding |
|-----------|--------|-------------|
| Replicate mouse study findings | Completed | Successfully validated the SEG-based approach on mouse ESC data |
| Process human SMART-Seq2 data | Completed | Generated count matrices from GSE75748 and E-MTAB-3929 |
| Implement refined SEG algorithm | Completed | Identified differentially stable genes between hESCs and developmental cells |
| Pathway enrichment analysis | Completed | Revealed enrichment in pluripotency pathways (OCT4, NANOG, SOX2 targets) |

---

## Part 1: Mouse Study Replication

### Dataset
- **E-MTAB-2600**: Mouse embryonic stem cells grown under three conditions (2i, a2i, serum)

### Quality Control Results

Cells were filtered using 3 MAD (median absolute deviation) thresholds on:
- Library size (total counts)
- Number of detected genes
- Mitochondrial gene proportion
- ERCC spike-in proportion

**QC Figures Generated:**

| Figure | Description |
|--------|-------------|
| `total_counts.png` | Distribution of library sizes across cells |
| `detected_genes.png` | Number of genes detected per cell |
| `mitochondrial_proportion.png` | Mitochondrial content percentage |
| `ERCC_proportion.png` | Spike-in proportion (technical control) |

### PCA Analysis

Using mouse stably expressed genes (mSEGs) from the scMerge package as features for PCA:

- Cells clustered distinctly by culture condition (2i, a2i, serum)
- SEG-based PCA showed clearer separation than using all genes
- Three clusters identified based on PC1 coordinates:
  - Cluster 0: PC1 between -6 and 4
  - Cluster 1: PC1 > 4
  - Cluster 2: PC1 < -6

**Key Finding:** The Seurat `RunPCA()` function is sensitive to feature selection, producing different results depending on whether features are explicitly or implicitly defined.

### Coefficient of Variation Analysis

Compared expression variability between:
- **SEGs**: Low coefficient of variation (CV), confirming stable expression
- **Non-SEGs**: Higher CV, demonstrating expected variability

This validated that the SEG list correctly identifies genes with stable expression patterns.

### NONO Gene Expression

The gene NONO (ENSMUSG00000031311) was examined across cell types:
- Differential expression observed between 2i, a2i, and serum conditions
- Demonstrates how a stably expressed gene can show condition-specific patterns

---

## Part 2: Human SMART-Seq2 Data Processing

### Datasets Processed

| Dataset | Source | Cell Type | Samples |
|---------|--------|-----------|---------|
| E-MTAB-3929 | ArrayExpress | Early human development (embryonic day 7) | Multiple developmental stages |
| GSE75748 | GEO | Undifferentiated human embryonic stem cells | hESC samples |

### Processing Pipeline

1. **Quality Control**: FastQC + MultiQC
2. **Trimming**: Trimmomatic (CROP:50)
3. **Alignment**: STAR to GRCh38 (Ensembl release 112)
4. **Quantification**: featureCounts
5. **Matrix Cleaning**: Column name standardization, metadata column removal

### Output Files
- `E-MTAB-3929_counts.csv`: Count matrix for developmental cells
- `GSE75748_counts.csv`: Count matrix for hESCs
- Corresponding `_targets.csv` metadata files

---

## Part 3: Human SEG Analysis Results

### Experimental Comparison

| Condition 1 | Condition 2 |
|-------------|-------------|
| E-MTAB-3929: Embryonic day 7 cells | GSE75748: Undifferentiated hESCs |
| Early developmental state | Pluripotent state |

### Quality Control

Both datasets underwent filtering:
- Low library size cells removed (3 MAD, lower)
- Low gene detection cells removed (3 MAD, lower)
- High mitochondrial content cells removed (3 MAD, higher)

### SEG Identification

Using the `scSEGIndex()` function from scMerge:

1. Calculated stability index (segIdx) for all genes in both datasets
2. Identified top 1000 most stable genes in each condition
3. Created union list of unique stable genes between conditions
4. Computed differential stability scores (segIdx_one - segIdx_two)

### Key Metrics

| Metric | E-MTAB-3929 | GSE75748 |
|--------|-------------|----------|
| Total genes | ~60,000 | ~60,000 |
| Protein-coding genes | Major subset | Major subset |
| Top 1000 SEGs | Identified | Identified |
| Unique stable genes | Determined from union list |

### Gene Set Enrichment Analysis (FGSEA)

Performed enrichment analysis using:
- **Database**: MSigDB v7.5.1 (Molecular Signatures Database)
- **Gene sets**: All pathways + Hallmark gene sets + KEGG pathways
- **Parameters**: minSize=5, maxSize=1000, nPermSimple=100,000

### Enrichment Results

#### Pathways Enriched in Developmental Cells (E-MTAB-3929)

Genes that are MORE stable in embryonic day 7 cells compared to hESCs showed enrichment in:
- Differentiation-related pathways
- Developmental signaling cascades

#### Pathways Enriched in hESCs (GSE75748)

Genes that are MORE stable in undifferentiated hESCs showed enrichment in:
- **BENPORATH_ES_1**: Embryonic stem cell signature genes
- **BENPORATH_NANOG_TARGETS**: NANOG transcription factor targets
- **BENPORATH_OCT4_TARGETS**: OCT4 (POU5F1) transcription factor targets
- **BENPORATH_SOX2_TARGETS**: SOX2 transcription factor targets

**Biological Interpretation**: The pluripotency master regulators (OCT4, NANOG, SOX2) maintain stable expression in undifferentiated hESCs, but their target genes become destabilized as cells differentiate.

### Specific Genes of Interest

| Gene | Ensembl ID | Role | Finding |
|------|------------|------|---------|
| NANOG | ENSG00000111704 | Pluripotency TF | Differential stability between conditions |
| SOX2 | ENSG00000181449 | Pluripotency TF | Key regulator with changed stability |
| PRDM14 | ENSG00000147596 | Pluripotency regulator | Tracked across conditions |
| POU5F1 (OCT4) | ENSG00000204531 | Master pluripotency TF | Core factor in enrichment |
| SFRP1 | ENSG00000104332 | WNT signaling | Differential expression observed |
| TCF7L1 | ENSG00000152284 | WNT pathway | Changed between states |

---

## Generated Output Files

### Analysis Tables

| File | Description |
|------|-------------|
| `original_EMTAB3929_SEG_table.csv` | Full SEG index for E-MTAB-3929 |
| `original_GSE75748_SEG_table.csv` | Full SEG index for GSE75748 |
| `unionlist_EMTAB3929_SEG_table.csv` | Filtered SEG table with differential scores |
| `unionlist_GSE75748_SEG_table.csv` | Filtered SEG table with differential scores |

### Enrichment Results

| File | Description |
|------|-------------|
| `EMTAB3929vsGSE75748_fgsea.csv` | Full GSEA results (developmental perspective) |
| `GSE75748vsEMTAB3929_fgsea.csv` | Full GSEA results (hESC perspective) |
| `EMTAB3929vsGSE75748_fgsea_hallmark.csv` | Hallmark pathway enrichment |
| `GSE75748vsEMTAB3929_fgsea_hallmark.csv` | Hallmark pathway enrichment |

### Figures

#### Quality Control Plots
- `figures/QC/total_counts.png`
- `figures/QC/detected_genes.png`
- `figures/QC/mitochondrial_proportion.png`
- `figures/QC/ERCC_proportion.png`

#### Analysis Plots
- `figures/other_plots/all_samples_pca.png` - PCA visualization of all cells
- `figures/other_plots/cluster_wise_pca.png` - PCA with cluster assignments
- `figures/other_plots/seg_cv.png` - Coefficient of variation comparison
- `figures/other_plots/nono_violin.png` - NONO expression across conditions

#### Enrichment Plots
- `combined_enrichment_plot.png` - OCT4, NANOG, SOX2 target enrichment

#### Heatmaps
- `E-MTAB-3929_Most_Stable.png` - Top 20 most stable genes (developmental)
- `E-MTAB-3929_Least_Stable.png` - Top 20 least stable genes (developmental)
- `GSE75748_Most_Stable.png` - Top 20 most stable genes (hESC)
- `GSE75748_Least_Stable.png` - Top 20 least stable genes (hESC)
- Gene-specific heatmaps: SFRP1, TCF7L1, POU5F1

---

## Conclusions

### Main Findings

1. **Method Validation**: The SEG-based approach successfully identified biologically meaningful differences between cell types, as demonstrated by replication of mouse study findings.

2. **Pluripotency Networks**: Human embryonic stem cells maintain stable expression of genes within the OCT4-NANOG-SOX2 regulatory network. These genes become destabilized during differentiation to embryonic day 7.

3. **Transcriptional Rewiring**: The transition from pluripotent to differentiated state involves fundamental changes in gene expression stability, not just expression levels.

4. **Pathway-Level Insights**: GSEA revealed that entire functional modules (pathways) show coordinated changes in stability, suggesting that the SEG-based approach captures higher-order biological organization.

### Technical Insights

1. **PCA Sensitivity**: Seurat's `RunPCA()` requires careful attention to feature selection for reproducible results.

2. **QC Importance**: MAD-based outlier detection effectively removes poor-quality cells while preserving biological variation.

3. **Cross-Dataset Comparison**: The union list approach enables meaningful comparison of stability indices across datasets with different gene detection rates.

### Biological Significance

The results support the model that:
- Pluripotent cells maintain a core transcriptional program with high stability
- Differentiation involves progressive destabilization of pluripotency-associated genes
- Gene stability (not just expression level) is a meaningful biological feature that reveals regulatory relationships

---

## Future Directions

1. Apply the algorithm to additional developmental timepoints
2. Integrate with chromatin accessibility data (ATAC-seq)
3. Validate key findings with independent datasets
4. Develop automated pipeline for broader application

---

*Document generated: 2026-02-04*
