# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Seurat")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")

library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)

# Load your Seurat object (assuming your Seurat object is already created and named)
# non_subsetted_all_genes_as_features <- readRDS("path_to_your_seurat_object.rds")

# Perform clustering (assuming PCA and ElbowPlot have been done and 10 PCs were selected)
non_subsetted_all_genes_as_features <- FindNeighbors(non_subsetted_all_genes_as_features, dims = 1:10)
non_subsetted_all_genes_as_features <- FindClusters(non_subsetted_all_genes_as_features, resolution = 0.5)

# Run UMAP for visualization
non_subsetted_all_genes_as_features <- RunUMAP(non_subsetted_all_genes_as_features, dims = 1:10)
DimPlot(non_subsetted_all_genes_as_features, reduction = "umap", label = TRUE) + ggtitle("UMAP Clustering")

# Run t-SNE for visualization (optional)
non_subsetted_all_genes_as_features <- RunTSNE(non_subsetted_all_genes_as_features, dims = 1:10)
DimPlot(non_subsetted_all_genes_as_features, reduction = "tsne", label = TRUE) + ggtitle("t-SNE Clustering")

# Identify cluster markers
cluster_markers <- FindAllMarkers(non_subsetted_all_genes_as_features, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, "cluster_markers.csv")

top10 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


# Annotate clusters manually (example)
non_subsetted_all_genes_as_features <- RenameIdents(non_subsetted_all_genes_as_features, `0` = "Cell Type A", `1` = "Cell Type B", `2` = "Cell Type C", `3` = "Cell Type D", `4` = "Cell Type E")

# Visualize marker expression
FeaturePlot(non_subsetted_all_genes_as_features, features = c("gene1", "gene2")) + ggtitle("Feature Plot")
VlnPlot(non_subsetted_all_genes_as_features, features = c("gene1", "gene2")) + ggtitle("Violin Plot")

# Perform differential expression analysis
cluster_0_vs_1 <- FindMarkers(non_subsetted_all_genes_as_features, ident.1 = 0, ident.2 = 1, min.pct = 0.25)
write.csv(cluster_0_vs_1, "cluster_0_vs_1.csv")

# Functional enrichment analysis using clusterProfiler
gene_list <- cluster_markers %>% filter(cluster == 0) %>% pull(gene)
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
go_results <- enrichGO(entrez_ids$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH")
dotplot(go_results) + ggtitle("GO Enrichment Analysis")

# Trajectory analysis (optional, using Monocle)
# library(monocle)
# cds <- as.cell_data_set(non_subsetted_all_genes_as_features)
# cds <- learn_graph(cds)
# plot_cells(cds, color_cells_by = "pseudotime") + ggtitle("Trajectory Analysis")

# Save the Seurat object and analysis results
saveRDS(non_subsetted_all_genes_as_features, file = "non_subsetted_all_genes_as_features.rds")