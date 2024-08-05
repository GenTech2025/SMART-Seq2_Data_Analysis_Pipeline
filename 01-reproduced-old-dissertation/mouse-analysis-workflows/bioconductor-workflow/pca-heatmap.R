library(Seurat)
library(dplyr)
library(pheatmap)

# Normalize columns using the provided function
normalize_columns <- function(x) {
  (x - mean(x)) / sd(x)
}

# Extract the PCA loading
subsetted_no_features_defined_loadings <- Loadings(subsetted_no_features_defined, reduction = "pca")[,1:2]
subsetted_no_features_defined_embeddings <- Embeddings(subsetted_no_features_defined, reduction = "pca")[,1:2]


# Convert loading to data frames
subsetted_no_features_defined_loadings <- as.data.frame(subsetted_no_features_defined_loadings)
subsetted_no_features_defined_embeddings <- as.data.frame(subsetted_no_features_defined_embeddings)


# Add a new column to the loading data frames which is the sum of first two principal components
subsetted_no_features_defined_loadings$SumOfPCs <- abs(subsetted_no_features_defined_loadings$PC_1) + abs(subsetted_no_features_defined_loadings$PC_2)
subsetted_mSEGs_as_features_loadings$SumOfPCs <- abs(subsetted_mSEGs_as_features_loadings$PC_1) + abs(subsetted_mSEGs_as_features_loadings$PC_2)


# Arrange the dataframes in descending order of SumOfPCs values
subsetted_no_features_defined_loadings <- subsetted_no_features_defined_loadings %>%
  arrange(desc(SumOfPCs))

subsetted_mSEGs_as_features_loadings <- subsetted_mSEGs_as_features_loadings %>%
  arrange(desc(SumOfPCs))

# Identify the top genes from the loadings dataframes based on SumOfPCs
top_genes_no_features <- rownames(subsetted_no_features_defined_loadings %>%
                                    arrange(desc(SumOfPCs)) %>%
                                    head(25))

top_genes_mSEGs <- rownames(subsetted_mSEGs_as_features_loadings %>%
                              arrange(desc(SumOfPCs)) %>%
                              head(25))

# Subset the Seurat objects based on the top genes
subsetted_no_features_defined_subset <- subset(subsetted_no_features_defined, features = top_genes_no_features)
subsetted_mSEGs_as_features_subset <- subset(subsetted_mSEGs_as_features, features = top_genes_mSEGs)

# Extract the log-normalized counts
lognorm_counts_no_features <- as.matrix(GetAssayData(subsetted_no_features_defined_subset, layer = "data"))
lognorm_counts_mSEGs <- as.matrix(GetAssayData(subsetted_mSEGs_as_features_subset, layer = "data"))

# Normalize columns using the provided function
lognorm_counts_no_features <- t(apply(lognorm_counts_no_features, 1, normalize_columns))
lognorm_counts_mSEGs <- t(apply(lognorm_counts_mSEGs, 1, normalize_columns))

# Transpose the matrices to have genes as columns and samples as rows
lognorm_counts_no_features <- t(lognorm_counts_no_features)
lognorm_counts_mSEGs <- t(lognorm_counts_mSEGs)

# Extract cell type annotations
cell_types_no_features <- as.factor(subsetted_no_features_defined_subset@meta.data$Type)
cell_types_mSEGs <- as.factor(subsetted_mSEGs_as_features_subset@meta.data$Type)

# Create heatmaps using pheatmap
# Heatmap for subsetted_no_features_defined_subset
# Define a custom color palette
custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Generate the heatmap with custom colors
pheatmap(lognorm_counts_no_features,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         annotation_row = data.frame(CellType = cell_types_no_features),
         annotation_colors = list(CellType = c("2i" = "purple", "a2i" = "orange", "serum" = "darkgreen")),
         main = "LogNormalized Counts No Features",
         show_rownames = FALSE,
         color = custom_colors)


# Define a custom color palette with dark red and dark blue
custom_colors <- colorRampPalette(c("darkblue", "white", "darkred"))(100)

# Generate the heatmap with custom colors and ensure the cell type annotation is shown
pheatmap(lognorm_counts_mSEGs,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_row = data.frame(CellType = cell_types_mSEGs),
         annotation_colors = list(CellType = c("2i" = "purple", "a2i" = "orange", "serum" = "darkgreen")),
         main = "LogNormalized Counts mSEGs",
         show_rownames = FALSE, # Do not show row names
         color = custom_colors)

# Heatmap of all genes as a reference
# Extract the normalized counts matrix
lognorm_counts <- as.matrix(GetAssayData(subsetted_no_features_defined, slot = "data"))

# Normalize columns using the provided function
lognorm_counts <- t(apply(lognorm_counts, 1, normalize_columns))

# Transpose the matrices to have genes as columns and samples as rows
lognorm_counts <- t(lognorm_counts)

# Check for NA/NaN/Inf values and handle them
lognorm_counts[is.na(lognorm_counts)] <- 0  # Replace NA values with 0
lognorm_counts[is.nan(lognorm_counts)] <- 0  # Replace NaN values with 0
lognorm_counts[is.infinite(lognorm_counts)] <- 0  # Replace Inf values with 0

# Generate the heatmap for all genes
pheatmap(lognorm_counts,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_row = data.frame(CellType = cell_types_no_features),
         annotation_colors = list(CellType = c("2i" = "purple", "a2i" = "orange", "serum" = "darkgreen")),
         main = "LogNormalized Counts No Features")
