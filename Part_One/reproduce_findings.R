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

# Subset the seurat object
so_norm_fil_mSEG <- subset(so_norm_fil, features=mSEG)

# Plot the PCA results
so_norm_fil <- RunPCA(so_norm_fil, features = rownames(so_norm_fil), verbose = FALSE)
so_norm_fil_mSEG <- RunPCA(so_norm_fil_mSEG, features = mSEG, verbose = FALSE)

# Plot and save the PCA results
p1 <- DimPlot(so_norm_fil, reduction = "pca", group.by = "Type") + ggtitle("All genes as features")
p2 <- DimPlot(so_norm_fil_mSEG, reduction = "pca", group.by = "Type") + ggtitle("Mouse Stably Expressed Genes as Features")
pca_compare <- cowplot::plot_grid(p1, p2, ncol = 2)

# Create subsets of the seurat object according to different cell types
so_raw_copy <- so_raw

# Normalize the data using seurat function
so_raw_copy <- NormalizeData(so_raw_copy, verbose = FALSE)
so_raw_copy <- ScaleData(so_raw_copy, verbose = FALSE)
so_raw_copy <- FindVariableFeatures(so_raw_copy, selection.method = "vst", verbose = FALSE)
variable_gene <- VariableFeatures(so_raw_copy)
so_seg_all <- subset(so_raw_copy, features = mSEG)

# Subset the seurat objects by cell type
so_raw_copy_serum <- subset(so_seg_all, Type == "serum")
so_raw_copy_2i <- subset(so_seg_all, Type == "2i")
so_raw_copy_a2i <- subset(so_seg_all, Type == "a2i")

# Compute CV for SEGs
so_seg_cv <- apply(as.matrix(GetAssayData(so_seg_all, layer = "data", assay = "RNA")),1,coefficient_of_variation)

# Determine remaining genes and reference genes
remaining_genes <- setdiff(rownames(so_raw_copy), mSEG)
reference_genes <- sample(remaining_genes, size = 722)

# Extract expression data for reference genes
reference_expr <- as.data.frame(GetAssayData(so_raw_copy, layer = "data", assay = "RNA"))
reference_expr <- reference_expr[reference_genes, ]

# Compute CV for nonSEGs
so_nonseg_cv <- apply(as.matrix(reference_expr), 1, coefficient_of_variation)

# Prepare and Plot the data
seg_plot <- data.frame(Gene = rownames(so_seg_all), CV = so_seg_cv)
nonseg_plot <- data.frame(Gene = rownames(reference_expr), CV = so_nonseg_cv)

cv_p1 <- ggplot(seg_plot, aes(x = Gene, y = CV)) +
  geom_point() +
  labs(title = "CV of SEG genes") +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  ylim(0, 30)

cv_p2 <- ggplot(nonseg_plot, aes(x = Gene, y = CV)) +
  geom_point() +
  labs(title = "CV of other genes") +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  ylim(0, 30)

cv_compare <- cowplot::plot_grid(cv_p1, cv_p2, ncol = 2)

# Save the plots
png(filename = "./figures/other_plots/seg_cv.png", width = 8, height = 8, units = "in", res = 300)

# Run PCA in the original seurat object which has all the cells
so_seg_all <- RunPCA(so_seg_all, assay = "RNA", verbose = FALSE)
pca_embedding <- as.data.frame(Embeddings(so_seg_all, reduction = "pca"))

png(filename = "./figures/other_plots/all_samples_pca.png", width = 8, height = 6, units = "in", res = 300)
DimPlot(so_seg_all, reduction = "pca", group.by = "Type") +
  ggtitle("PCA for all cells (feature = SEG)")
dev.off()

# Determine which cluster each cell belongs to
so_seg_all$seurat_clusters <- ifelse(
  pca_embedding$PC_1 > -6 & pca_embedding$PC_1 <=4, "cluster0", ifelse(pca_embedding$PC_1 > 4, "cluster1", "cluster2")
)

# set active identity class to new clusters
Idents(so_seg_all) <- factor(so_seg_all$seurat_clusters)

# PLot clusters for all cells with vertical lines
png("./figures/other_plots/cluster_wise_pca.png", width = 8, height = 6, units = "in", res = 300)
DimPlot(so_seg_all, reduction = "pca", group.by = "seurat_clusters") +
  ggtitle("Cluster for all cells") +
  geom_vline(xintercept = 4, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = -6, linetype = "dashed", color = "blue")
dev.off()

# Create a violin plot to visualize the differential expression values of Nono across different cell types
# Ensure that only the specific cell types are included
selected_cell_types <- subset(so_raw_copy, subset = Type %in% c("2i", "a2i", "serum"))

# Extract the expression data for the gene of interest and the cell type information
gene_expression_data <- FetchData(selected_cell_types, vars = c("ENSMUSG00000031311", "Type"))

# Filter out any unwanted cell types explicitly (if necessary)
gene_expression_data <- gene_expression_data[gene_expression_data$Type %in% c("2i", "a2i", "serum"), ]

# Create the enhanced box plot (Note that the original code of violin plot could not be found)
nono_violin <- ggplot(gene_expression_data, aes(x = Type, y = ENSMUSG00000031311, fill = Type)) +
                        geom_violin(trim = TRUE, alpha = 0.6, color = "black") +  # Trim for symmetry
                        stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red", fill = "white") +  # Add mean points
                        scale_fill_manual(values = c("2i" = "red", "a2i" = "blue", "serum" = "green")) +  # Updated colors
                        labs(
                            title = "Expression of Nono Across Cell Types",
                            x = "Cell Type",
                            y = "Expression Level",
                            fill = "Cell Type"
                        ) +
                        theme_minimal(base_size = 15) +
                        theme(
                            plot.title = element_text(hjust = 0.5, face = "bold"),
                            axis.title.x = element_text(face = "bold"),
                            axis.title.y = element_text(face = "bold"),
                            legend.position = "right",
                            panel.grid.major = element_line(color = "grey80"),
                            axis.text.x = element_text(angle = 45, hjust = 1)
                        )

# Print the enhanced plot
png(filename = "./figures/other_plots/nono_violin.png", width = 8, height = 6, units = "in", res = 300)
print(nono_violin)
dev.off()
