###### THE PCA OBJECT SHOULD EXIST IN THE ENVIRONMENT###########


# Extract embeddings and loadings for both datasets
embeddings_1 <- Embeddings(subsetted_no_features_defined, reduction = "pca")
loadings_1 <- as.data.frame(subsetted_no_features_defined@reductions$pca@feature.loadings)

embeddings_2 <- Embeddings(subsetted_mSEGs_as_features, reduction = "pca")
loadings_2 <- as.data.frame(subsetted_mSEGs_as_features@reductions$pca@feature.loadings)


# Calculate combined contribution for PC1 and PC2
loadings_1$combined <- abs(loadings_1$PC_1) + abs(loadings_1$PC_2)
loadings_2$combined <- abs(loadings_2$PC_1) + abs(loadings_2$PC_2)

# Get top 10 genes based on the combined contribution
top_genes_1 <- rownames(loadings_1[order(loadings_1$combined, decreasing = TRUE), ])[1:10]
top_genes_2 <- rownames(loadings_2[order(loadings_2$combined, decreasing = TRUE), ])[1:10]

# Combine top genes and ensure there are 20 unique genes
all_genes <- unique(c(top_genes_1, top_genes_2))

# If there are fewer than 20 genes (due to overlaps), add additional genes to make 20
if (length(all_genes) < 20) {
  additional_genes_1 <- setdiff(rownames(loadings_1[order(loadings_1$combined, decreasing = TRUE), ]), all_genes)
  additional_genes_2 <- setdiff(rownames(loadings_2[order(loadings_2$combined, decreasing = TRUE), ]), all_genes)
  all_genes <- unique(c(all_genes, additional_genes_1, additional_genes_2))
  all_genes <- all_genes[1:20]
}

# Map each gene to its origin PCA plot
gene_origin <- sapply(all_genes, function(gene) {
  origin <- c()
  if (gene %in% top_genes_1) origin <- c(origin, "PCA Plot 1")
  if (gene %in% top_genes_2) origin <- c(origin, "PCA Plot 2")
  paste(origin, collapse = ", ")
})

top_genes <- data.frame(Gene = all_genes, Origin = gene_origin, stringsAsFactors = FALSE)


# Align matrices to have the same set of genes
valid_genes_1 <- intersect(top_genes$Gene, rownames(loadings_1))
valid_genes_2 <- intersect(top_genes$Gene, rownames(loadings_2))

# Recreate top_genes to ensure all genes from both PCA plots are included
valid_genes <- unique(c(valid_genes_1, valid_genes_2))

# Calculate contribution matrices for both PCA plots using the combined list of valid genes
contribution_matrix_1 <- embeddings_1[, c("PC_1", "PC_2")] %*% t(loadings_1[valid_genes, c("PC_1", "PC_2")])
contribution_matrix_2 <- embeddings_2[, c("PC_1", "PC_2")] %*% t(loadings_2[valid_genes, c("PC_1", "PC_2")])


# Function to normalize each column
normalize_columns <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Apply normalization
normalized_matrix_1 <- apply(contribution_matrix_1, 2, normalize_columns)
normalized_matrix_2 <- apply(contribution_matrix_2, 2, normalize_columns)


# Ensure both normalized matrices have the same rows (samples)
rownames(normalized_matrix_1) <- rownames(embeddings_1)
rownames(normalized_matrix_2) <- rownames(embeddings_2)

# Combine the normalized matrices into a single matrix by columns
combined_matrix <- cbind(normalized_matrix_1, normalized_matrix_2)

# Ensure the row names (sample identifiers) are unique
rownames(combined_matrix) <- rownames(normalized_matrix_1)

# Verify lengths
print(paste("Number of rows in combined_matrix:", nrow(combined_matrix)))
print(paste("Number of cell types:", length(c(subsetted_no_features_defined@meta.data$Type, subsetted_mSEGs_as_features@meta.data$Type))))
print(paste("Number of columns in combined_matrix:", ncol(combined_matrix)))
print(paste("Number of top genes:", nrow(top_genes)))

# Create annotation data frame for cell types
cell_types <- subsetted_no_features_defined@meta.data$Type

# Ensure the lengths match
if (nrow(combined_matrix) != length(cell_types)) {
  stop("Mismatch between number of rows in combined_matrix and length of cell_types")
}

annotation_df <- data.frame(CellType = cell_types, stringsAsFactors = FALSE)
rownames(annotation_df) <- rownames(combined_matrix)

# Generate heatmap
library(ComplexHeatmap)

# Update top_genes to ensure all genes are included and correctly ordered
top_genes_combined <- data.frame(
  Gene = colnames(combined_matrix),
  Origin = sapply(colnames(combined_matrix), function(g) {
    if (g %in% top_genes$Gene) {
      return(top_genes$Origin[match(g, top_genes$Gene)])
    } else {
      return("Both PCA Plots")
    }
  }),
  stringsAsFactors = FALSE
)

png("./combined_heatmap.png", width = 12, height = 10, units = "in", res = 500)
Heatmap(combined_matrix, name = "Normalized Contribution",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = TRUE, show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "Contribution"),
        top_annotation = HeatmapAnnotation(
          PlotOrigin = anno_simple(top_genes_combined$Origin, col = c("PCA Plot 1" = "blue", "PCA Plot 2" = "red", "Both PCA Plots" = "purple"))
        ),
        right_annotation = rowAnnotation(CellType = annotation_df$CellType,
                                         col = list(CellType = c("a2i" = "orange", "2i" = "purple", "serum" = "green"))))
dev.off()