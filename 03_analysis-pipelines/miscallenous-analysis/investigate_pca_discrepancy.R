# Data directory
data_dir <- "/home/s2599932/study-datasets/mouse-dataset/E-MTAB-2600/single_es_tximport.rds"
annotation_dir <- "/home/s2599932/study-datasets/mouse-dataset/E-MTAB-2600/E-MTAB-2600.targets.txt"

# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(scran)
library(scater)
library(biomaRt)
library(scMerge)
library(cowplot)

# Load data
raw <- readRDS(data_dir)
annotation <- read.table(annotation_dir, sep = "\t", header = TRUE, row.names = "ERR")
annotation <- annotation[match(colnames(raw$counts), rownames(annotation)), ]

data("segList_ensemblGeneID", package = "scMerge")
mSEG <- segList_ensemblGeneID$mouse$mouse_scSEG

# Create SCE object
sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(raw$counts)),
  colData = annotation
)

# Perform Quality Control
is_spike <- grepl("^ERCC", rownames(sce))
mitogenes <- "ENSMUSG00000064336|ENSMUSG00000064337|ENSMUSG00000064338|ENSMUSG00000064339|ENSMUSG00000064340|ENSMUSG00000064341|ENSMUSG00000064342|ENSMUSG00000064343|ENSMUSG00000064344|ENSMUSG00000064345|ENSMUSG00000064346|ENSMUSG00000064347|ENSMUSG00000064348|ENSMUSG00000064349|ENSMUSG00000064350|ENSMUSG00000064351|ENSMUSG00000064352|ENSMUSG00000064353|ENSMUSG00000064354|ENSMUSG00000064355|ENSMUSG00000064356|ENSMUSG00000064357|ENSMUSG00000064358|ENSMUSG00000064359|ENSMUSG00000064360|ENSMUSG00000064361|ENSMUSG00000065947|ENSMUSG00000064363|ENSMUSG00000064364|ENSMUSG00000064365|ENSMUSG00000064366|ENSMUSG00000064367|ENSMUSG00000064368|ENSMUSG00000064369|ENSMUSG00000064370|ENSMUSG00000064371|ENSMUSG00000064372"
is_mito <- grepl(mitogenes, rownames(sce))
sce$qc_stats <- perCellQCMetrics(sce, subsets=list(ERCC=is_spike, Mt=is_mito))

# Calculate the QC metrics based on which cells will be retained or eliminated
libsize_drop <- isOutlier(sce$qc_stats$sum, nmads = 3, type = "lower", log = TRUE)
feature_drop <- isOutlier(sce$qc_stats$detected, nmads = 3, type = "lower", log = TRUE)
mito_drop <- isOutlier(sce$qc_stats$subsets_Mt_percent, nmads = 3, type = "higher")
spike_drop <- isOutlier(sce$qc_stats$subsets_ERCC_percent, nmads = 3, type = "higher")

# Save a copy of the unfiltered SCE object
sce_copy <- sce

# Drop the samples from the SCE object that do not meet the QC standards
sce_fil <- sce[, !(libsize_drop | feature_drop | mito_drop | spike_drop)]
filtered_cells <- colnames(sce_fil)

# Normalize the expression data using deconvolution normalization
clusters <- quickCluster(sce_fil)
sce_fil <- computeSumFactors(sce_fil, clusters = clusters)
sce_fil <- logNormCounts(sce_fil)

# Run and Plot the PCA results using RunPCA()
# Create Seurat Object
log_counts <- assay(sce_fil, "logcounts")

# Creating the Seurat object
so <- CreateSeuratObject(counts = log_counts)
so <- AddMetaData(so, metadata = annotation)

DefaultAssay(so) <- "RNA"
so <- SetAssayData(so, layer = "data", new.data = log_counts)

so_fil <- so[, filtered_cells]

# Only include single cell samples
# Convert the metadata to a tibble for filtering
metadata <- so_fil@meta.data %>%
  rownames_to_column(var = "cell")

# Filter the metadata to only include the desired cell types
filtered_metadata <- metadata %>%
  filter(Type %in% c("2i", "a2i", "serum"))

# Subset the Seurat object using the filtered metadata
so_fil <- subset(so_fil, cells = filtered_metadata$cell)

# Scale the data
so_fil <- ScaleData(so_fil, verbose = FALSE, features = rownames(so_fil))
so_fil <- FindVariableFeatures(so_fil, selection.method = "vst", verbose = FALSE)

# Make copies of the seurat objects
so_non_subsetted_no_features_defined <- so_fil
so_non_subsetted_all_genes_as_features <- so_fil
#so_non_subsetted_mSEGs_as_features <- so_fil
so_subsetted_no_features_defined <- so_fil[mSEG,]
so_subsetted_mSEGs_as_features <- so_fil[mSEG,]

# Define the PCA plotting function
plot_pca_seurat <- function(seurat_obj, features, metadata_column, title) {
  # Run PCA on the Seurat object using the provided features
  seurat_obj <- RunPCA(seurat_obj, features = features, verbose = FALSE)
  
  # Extract PCA embeddings and metadata, and combine them
  pca_df <- Embeddings(seurat_obj, reduction = "pca")[, 1:2] %>%
    as.data.frame() %>%
    rownames_to_column(var = "cell") %>%
    left_join(seurat_obj@meta.data %>%
                rownames_to_column(var = "cell"),
              by = "cell")
  
  # Calculate the variance explained by the first two principal components (for labels)
  explained_variance <- Stdev(seurat_obj, reduction = "pca")^2 / sum(Stdev(seurat_obj, reduction = "pca")^2) * 100
  pc1 <- round(explained_variance[1], 2)
  pc2 <- round(explained_variance[2], 2)
  
  # Plot with ggplot2 using the specified parameters
  pca_plot <- ggplot(pca_df, aes_string(x = "PC_1", y = "PC_2", color = metadata_column)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_brewer(palette = "Set1") +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold", hjust = 0.5, size = 12),
      legend.position = "right",
      legend.key.size = unit(0.8, "cm"),
      legend.text = element_text(size = 10),
      plot.margin = margin(20, 20, 20, 20)  # Add margins (top, right, bottom, left)
    ) +
    labs(
      title = title,
      x = paste("PC_1 (", pc1, "%)", sep = ""),
      y = paste("PC_2 (", pc2, "%)", sep = ""),
      color = "Cell Type"
    )
  
  return(pca_plot)
}

# Plot 1
p1 <- plot_pca_seurat(seurat_obj = so_non_subsetted_no_features_defined, features = NULL, metadata_column = "Type", title = "RunPCA() + Non Subsetted + No features defined")
print(p1)

# Plot 2
p2 <- plot_pca_seurat(seurat_obj = so_non_subsetted_all_genes_as_features, features = rownames(so_non_subsetted_all_genes_as_features), metadata_column = "Type", title = "RunPCA() + Non subsetted + All genes as features")
print(p2)

# Plot 3
p3 <- plot_pca_seurat(seurat_obj = so_subsetted_no_features_defined, features = NULL, metadata_column = "Type", title = "RunPCA() + Subsetted + No features defined")
print(p3)

# Plot 4
p4 <- plot_pca_seurat(seurat_obj = so_subsetted_mSEGs_as_features, features = mSEG, metadata_column = "Type", title = "RunPCA() + Subsetted + mSEGs as features")
print(p4)

# Run and Plot the PCA results using runPCA()
# Only include the single cells
# Remove the bulk samples and only include single cell samples
metadata <- as.data.frame(colData(sce_fil)) %>%
  rownames_to_column(var = "cell")

filtered_metadata <- metadata %>%
  filter(Type %in% c("2i", "a2i", "serum"))

sce_fil <- sce_fil[, colnames(sce_fil) %in% filtered_metadata$cell]

mSEG_index <- intersect(mSEG, rownames(sce_non_subsetted_mSEGs_as_features))

# Copy the SCE objects
sce_non_subsetted_no_features_defined <- sce_fil
sce_non_subsetted_all_genes_as_features <- sce_fil
sce_subsetted_no_features_defined <- sce_fil[mSEG_index,]
sce_subsetted_mSEGs_as_features <- sce_fil[mSEG_index,]

# Define the PCA plotting function
plot_pca_scater <- function(sce_obj, exprs_values = "logcounts", subset_row = NULL, scale = TRUE, metadata_column, title) {
  # Run PCA on the SingleCellExperiment object
  sce_obj <- runPCA(sce_obj, exprs_values = exprs_values, subset_row = subset_row, scale = scale)
  
  # Extract PCA embedding
  pca_embeddings <- reducedDim(sce_obj, "PCA")[, 1:2]
  pca_df <- as.data.frame(pca_embeddings)
  colnames(pca_df) <- c("PC_1", "PC_2")
  pca_df <- rownames_to_column(pca_df, var = "cell")
  
  # Combine with metadata
  metadata <- as.data.frame(colData(sce_obj))
  metadata <- rownames_to_column(metadata, var = "cell")
  pca_df <- left_join(pca_df, metadata, by = "cell")
  
  # Calculate the variance explained by the first two principal components
  explained_variance <- attr(reducedDim(sce_obj, "PCA"), "percentVar")
  pc1 <- round(explained_variance[1], 2)
  pc2 <- round(explained_variance[2], 2)
  
  # Plot with ggplot2
  pca_plot <- ggplot(pca_df, aes_string(x = "PC_1", y = "PC_2", color = metadata_column)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_brewer(palette = "Set1") +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold", hjust = 0.5, size = 12),
      legend.position = "right",
      legend.key.size = unit(0.6, "cm"),
      legend.text = element_text(size = 10),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    labs(
      title = title,
      x = paste("PC_1 (", pc1, "%)", sep = ""),
      y = paste("PC_2 (", pc2, "%)", sep = ""),
      color = "Cell Type"
    )
  
  return(pca_plot)
}

# Plot 1
scater_p1 <- plot_pca_scater(sce_obj =sce_non_subsetted_no_features_defined, exprs_values = "logcounts", subset_row = NULL, scale = TRUE, metadata_column = "Type", title = "runPCA() + Non substetted + No features defined" )
print(scater_p1)

# Plot 2
scater_p2 <- plot_pca_scater(sce_obj =sce_non_subsetted_all_genes_as_features, exprs_values = "logcounts", subset_row = rownames(sce_non_subsetted_all_genes_as_features), scale = TRUE, metadata_column = "Type", title = "runPCA() + Non substetted + All genes as features" )
print(scater_p2)

# Plot 3
scater_p3 <- plot_pca_scater(sce_obj =sce_subsetted_no_features_defined, exprs_values = "logcounts", subset_row = NULL, scale = TRUE, metadata_column = "Type", title = "runPCA() + Substetted + No features defined" )
print(scater_p3)

# Plot 4
scater_p4 <- plot_pca_scater(sce_obj =sce_subsetted_mSEGs_as_features, exprs_values = "logcounts", subset_row = mSEG_index, scale = TRUE, metadata_column = "Type", title = "runPCA() + Substetted + mSEGs as features" )
print(scater_p4)

# Run and Plot the PCA results using fixedPCA()
# Only include the single cells
# Remove the bulk samples and only include single cell samples
metadata <- as.data.frame(colData(sce_fil)) %>%
  rownames_to_column(var = "cell")

filtered_metadata <- metadata %>%
  filter(Type %in% c("2i", "a2i", "serum"))

sce_fil <- sce_fil[, colnames(sce_fil) %in% filtered_metadata$cell]

mSEG_index <- intersect(mSEG, rownames(sce_fil))

# Copy the SCE objects
scran_sce_non_subsetted_no_features_defined <- sce_fil
scran_sce_non_subsetted_all_genes_as_features <- sce_fil
scran_sce_subsetted_no_features_defined <- sce_fil[mSEG_index,]
scran_sce_subsetted_mSEGs_as_features <- sce_fil[mSEG_index,]

# Define the PCA plotting function
plot_pca_scran <- function(sce_obj, subset_row = NULL, assay_type = "logcounts", metadata_column, title) {
  # Perform PCA on the SingleCellExperiment object using fixedPCA
  sce_obj <- fixedPCA(sce_obj, subset.row = subset_row, assay.type = assay_type)
  
  # Extract PCA embedding
  pca_embeddings <- reducedDim(sce_obj, "PCA")[, 1:2]
  pca_df <- as.data.frame(pca_embeddings)
  colnames(pca_df) <- c("PC_1", "PC_2")
  pca_df <- rownames_to_column(pca_df, var = "cell")
  
  # Combine with metadata
  metadata <- as.data.frame(colData(sce_obj))
  metadata <- rownames_to_column(metadata, var = "cell")
  pca_df <- left_join(pca_df, metadata, by = "cell")
  
  # Calculate the variance explained by the first two principal components
  explained_variance <- attr(reducedDim(sce_obj, "PCA"), "percentVar")
  pc1 <- round(explained_variance[1], 2)
  pc2 <- round(explained_variance[2], 2)
  
  # Plot with ggplot2
  pca_plot <- ggplot(pca_df, aes_string(x = "PC_1", y = "PC_2", color = metadata_column)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_brewer(palette = "Set1") +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold", hjust = 0.5, size = 12),
      legend.position = "right",
      legend.key.size = unit(0.6, "cm"),
      legend.text = element_text(size = 10),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    labs(
      title = title,
      x = paste("PC_1 (", pc1, "%)", sep = ""),
      y = paste("PC_2 (", pc2, "%)", sep = ""),
      color = "Cell Type"
    )
  
  return(pca_plot)
}

# Plot 1
scran_p1 <- plot_pca_scran(sce_obj = scran_sce_non_subsetted_no_features_defined, subset_row = NULL, assay_type = "logcounts", metadata_column = "Type", title = "fixedPCA() + Non subsetted + No features defined")
print(scran_p1)

# Plot 2
scran_p2 <- plot_pca_scran(sce_obj = scran_sce_non_subsetted_all_genes_as_features, subset_row = rownames(scran_sce_non_subsetted_all_genes_as_features), assay_type = "logcounts", metadata_column = "Type", title = "fixedPCA() + Non subsetted + All genes as features")
print(scran_p2)

# Plot 3
scran_p3 <- plot_pca_scran(sce_obj = scran_sce_subsetted_no_features_defined, subset_row = NULL, assay_type = "logcounts", metadata_column = "Type", title = "fixedPCA() + Subsetted + No features defined")
print(scran_p3)

# Plot 4
scran_p4 <- plot_pca_scran(sce_obj = scran_sce_subsetted_mSEGs_as_features, subset_row = mSEG_index, assay_type = "logcounts", metadata_column = "Type", title = "fixedPCA() + Subsetted + mSEGs as features")
print(scran_p4)

# Create panel plots to facilitate the comparison between the three functions
# Panel plot for all genes as features or no features
all_genes_as_features_or_no_features <- plot_grid(p1, p2, scater_p1, scater_p2, scran_p1, scran_p2, nrow = 3, ncol = 2, labels = c("A", "B", "C", "D", "E", "F"))
print(all_genes_as_features_or_no_features)

# Panel plot for mSEGs as features or no features
mSEG_as_features_or_no_features <- plot_grid(p3, p4, scater_p3, scater_p4, scran_p3, scran_p4, nrow = 3, ncol = 2, labels = c("A", "B", "C", "D", "E", "F"))
print(mSEG_as_features_or_no_features)

# Define function for saving plots
save_high_quality_png <- function(plot, filename, width = 14, height = 14, dpi = 400) {
  ggsave(filename = filename, plot = plot, width = width, height = height, dpi = dpi, device = "png")
}

# Save the plots
save_high_quality_png(all_genes_as_features_or_no_features, filename = "./all_genes_PCA_panel_plot.png")
save_high_quality_png(mSEG_as_features_or_no_features, filename = "./mSEGs_PCA_panel_plot.png")

# Save the environment information
save.image(file = "mouse-pca.Rdata")
