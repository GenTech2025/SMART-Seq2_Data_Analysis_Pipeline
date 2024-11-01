# Reproduced analysis for the last years dissertation with some improvements and modifications

# Load packages
library(Seurat)
library(scMerge)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(corrplot)
library(scater)
library(scran)

# Define functions
coefficient_of_variation <- function(x){
  cv <- sd(x)/mean(x)
  return (cv)
}

# Load in data and annotation file
es_raw <- readRDS("/home/s2599932/mouse-dataset/E-MTAB-2600/single_es_tximport.rds")
es_targets <- read.table("/home/s2599932/mouse-dataset/E-MTAB-2600/E-MTAB-2600.targets.txt", sep = "\t", header = TRUE, row.names = "ERR")

# Store mSEGs defined in scMerge in a variable
data("segList_ensemblGeneID", package = "scMerge")
mSEG <- segList_ensemblGeneID$mouse$mouse_scSEG

# Create a SingleCellExperiment Object
es_sce <- SingleCellExperiment(assays=list(counts = es_raw$counts, abundance = es_raw$abundance))

# Identify mouse mitochondrial genes
is_spike <- grepl("^ERCC",rownames(es_sce))
mitogenes <- "ENSMUSG00000064336|ENSMUSG00000064337|ENSMUSG00000064338|ENSMUSG00000064339|ENSMUSG00000064340|ENSMUSG00000064341|ENSMUSG00000064342|ENSMUSG00000064343|ENSMUSG00000064344|ENSMUSG00000064345|ENSMUSG00000064346|ENSMUSG00000064347|ENSMUSG00000064348|ENSMUSG00000064349|ENSMUSG00000064350|ENSMUSG00000064351|ENSMUSG00000064352|ENSMUSG00000064353|ENSMUSG00000064354|ENSMUSG00000064355|ENSMUSG00000064356|ENSMUSG00000064357|ENSMUSG00000064358|ENSMUSG00000064359|ENSMUSG00000064360|ENSMUSG00000064361|ENSMUSG00000065947|ENSMUSG00000064363|ENSMUSG00000064364|ENSMUSG00000064365|ENSMUSG00000064366|ENSMUSG00000064367|ENSMUSG00000064368|ENSMUSG00000064369|ENSMUSG00000064370|ENSMUSG00000064371|ENSMUSG00000064372"
is_mito <- grepl(mitogenes, rownames(es_sce))

# Calculate quality control metrics
es_sce$qc_stats <- perCellQCMetrics(es_sce, subsets=list(ERCC=is_spike, Mt=is_mito))

# Visualize the different quality control plots
# Library sizes in millions
tiff("/home/s2599932/dissertation-main/plots/reproduced-old-dissertation/quality-control/library_sizes.tiff", width=7, height=7, units="in", res=300)
hist(es_sce$qc_stats$sum/1e6, 
     xlab="Library sizes (millions)", 
     main="Histogram of Library Sizes", 
     breaks=20, 
     col="grey80", 
     ylab="Number of cells")
dev.off()

# Number of detected genes
tiff("/home/s2599932/dissertation-main/plots/reproduced-old-dissertation/quality-control/detected_genes.tiff", width=7, height=7, units="in", res=300)
hist(es_sce$qc_stats$detected, 
     xlab="Number of detected genes", 
     main="Histogram of Detected Genes", 
     breaks=20, 
     col="grey80", 
     ylab="Number of cells")
dev.off()

# Mitochondrial gene proportion (%)
tiff("/home/s2599932/dissertation-main/plots/reproduced-old-dissertation/quality-control/mitochondial_proportion.tiff", width=7, height=7, units="in", res=300)
hist(es_sce$qc_stats$subsets_Mt_percent,
     xlab="Mitochondrial proportion (%)", 
     main="Histogram of Mitochondrial Proportion", 
     breaks=20, 
     col="grey80", 
     ylab="Number of cells")
dev.off()

# ERCC proportion (%)
tiff("/home/s2599932/dissertation-main/plots/reproduced-old-dissertation/quality-control/ERCC_proportion.tiff", width=7, height=7, units="in", res=300)
hist(es_sce$qc_stats$subsets_ERCC_percent,
     xlab="ERCC proportion (%)", 
     main="Histogram of ERCC Proportion", 
     breaks=20, 
     col="grey80", 
     ylab="Number of cells")
dev.off()

# Calculating the filtering metrics
libsize_drop <- isOutlier(es_sce$qc_stats$sum, nmads = 3, type = "lower", log = TRUE)
feature_drop <- isOutlier(es_sce$qc_stats$detected, nmads = 3, type = "lower", log = TRUE)
mito_drop <- isOutlier(es_sce$qc_stats$subsets_Mt_percent, nmads = 3, type = "higher")
spike_drop <- isOutlier(es_sce$qc_stats$subsets_ERCC_percent, nmads = 3, type = "higher")

# Remove low quality samples/cells from the SingleCellExperiment(SCE) object
es_sce_fil <- es_sce[,!(libsize_drop|feature_drop|mito_drop|spike_drop)]

# store the names of the samples/cells that passed QC
filtered_cells <- colnames(es_sce_fil)

