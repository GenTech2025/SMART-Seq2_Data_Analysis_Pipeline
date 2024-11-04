# Load libraries
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
  # Function to calculate coefficient of variation
  cv <- sd(x)/mean(x)
  return (cv)
}

save_histogram <- function(data, filename, path, xlab, title){
  # Function to save histogram
  png(filename = paste0(path, filename), width = 6, height = 6, units = "in", res = 300)
  hist(data, 
     xlab=xlab, 
     main=title, 
     breaks=20,
     col="grey80",
     ylab="Number of cells")
  dev.off()
}

# Load in the data and annotation file
es_targets <- read.table("./data/mouse-dataset/E-MTAB-2600/E-MTAB-2600.targets.txt", sep = "\t", header = TRUE, row.names = "ERR")
es_raw <- readRDS("./data/mouse-dataset/E-MTAB-2600/single_es_tximport.rds")

# Create a SingleCellExperiment object
es_sce <- SingleCellExperiment(assays=list(counts = es_raw$counts, abundance = es_raw$abundance))

# Load in the mouse stably expressed genes
data("segList_ensemblGeneID", package = "scMerge")
mSEG <- segList_ensemblGeneID$mouse$mouse_scSEG

# Connect to the Ensembl BioMart database
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Retrieve mitochondrial genes
mito_genes <- getBM(attributes = c("ensembl_gene_id"), 
                    filters = "chromosome_name", 
                    values = "MT", 
                    mart = ensembl)

# Create logical vectors for spike-in and mitochondrial genes
is_spike <- grepl("^ERCC", rownames(es_sce))
is_mito <- rownames(es_sce) %in% mito_genes$ensembl_gene_id

# Calculate the Quality Control Metrics
es_sce$qc_stats <- perCellQCMetrics(es_sce, subsets=list(ERCC=is_spike, Mt=is_mito))

# Save the QC histograms
save_histogram(es_sce$qc_stats$sum/1e6, "total_counts.png", "./figures/QC/", "Library sizes (millions)", "Histogram of Library Sizes")
save_histogram(es_sce$qc_stats$detected, "detected_genes.png", "./figures/QC/", xlab="Number of genes detected", "Histogram of Detected Genes")
save_histogram(es_sce$qc_stats$subsets_Mt_percent, "mitochondrial_proportion.png", "./figures/QC/", xlab="Mitochondrial Proportion (%)", "Histogram of Mitochondrial Proportion")
save_histogram(es_sce$qc_stats$subsets_ERCC_percent, "ERCC_proportion.png", "./figures/QC/", xlab="ERCC Proportion (%)", "Histogram of ERCC Proportion")

# Filter cells based on QC metrics
libsize_drop <- isOutlier(es_sce$qc_stats$sum, nmads = 3, type = "lower", log = TRUE)
feature_drop <- isOutlier(es_sce$qc_stats$detected, nmads = 3, type = "lower", log = TRUE)
mito_drop <- isOutlier(es_sce$qc_stats$subsets_Mt_percent, nmads = 3, type = "higher")
spike_drop <- isOutlier(es_sce$qc_stats$subsets_ERCC_percent, nmads = 3, type = "higher")
 
es_sce_copy <- es_sce # Save a non filtered copy

# Filtering the SingleCellExperiment(SCE) object
es_sce_fil <- es_sce[,!(libsize_drop|feature_drop|mito_drop|spike_drop)]
es_sce_fil_copy <- es_sce_fil # Save a filtered copy

# Save the filtered cells
filtered_cells <- colnames(es_sce_fil)

# Create a Seurat Object from the counts data
es_counts <- as.data.frame(es_raw$counts)
so_raw <- CreateSeuratObject(counts = es_counts)
so_raw <- AddMetaData(so_raw, metadata = es_targets)

# Filter the Seurat object to only include the cells/samples that passed QC
so_raw_fil <- so_raw[,filtered_cells]

# Normalize, select variable features and scale the data
so_norm_fil <- NormalizeData(so_raw_fil, verbose = FALSE)
so_norm_fil <- FindVariableFeatures(so_norm_fil, selection.method = "vst", verbose = FALSE)
so_norm_fil <- ScaleData(so_norm_fil,features = rownames(so_norm_fil), verbose = FALSE)

# Create copies of the seurat object
so_norm_fil_mSEG <- so_norm_fil

# Run PCA using all genes as features and using mouse stably expressed genes as features
so_norm_fil <- RunPCA(so_norm_fil, features = rownames(so_norm_fil), verbose = FALSE)
# Subset the seurat object
so_norm_fil_mSEG <- subset(so_norm_fil_mSEG, features=mSEG)
so_norm_fil_mSEG <- RunPCA(so_norm_fil_mSEG, features = mSEG, verbose = FALSE)

# Plot the results
p1<-DimPlot(so_norm_fil,reduction= "pca",group.by = "Type")+ ggtitle("All genes as features")
p2<-DimPlot(so_norm_fil_mSEG, reduction = "pca",group.by = "Type")+ ggtitle("Mouse stably expressed genes as features")
pca_compare <- cowplot::plot_grid(p1, p2, ncol = 2)