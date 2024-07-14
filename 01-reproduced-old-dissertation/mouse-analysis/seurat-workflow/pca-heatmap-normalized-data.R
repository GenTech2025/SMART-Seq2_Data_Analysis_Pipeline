library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# Ensure the relevant columns are numeric
subsetted_no_features_defined_loadings <- subsetted_no_features_defined_loadings %>%
  mutate(PC_1 = as.numeric(PC_1), PC_2 = as.numeric(PC_2), SumOfPCs = as.numeric(SumOfPCs))

subsetted_mSEGs_as_features_loadings <- subsetted_mSEGs_as_features_loadings %>%
  mutate(PC_1 = as.numeric(PC_1), PC_2 = as.numeric(PC_2), SumOfPCs = as.numeric(SumOfPCs))

# Identify the top 25 genes from the loadings dataframes based on SumOfPCs
top_genes_no_features <- subsetted_no_features_defined_loadings %>%
  arrange(desc(SumOfPCs)) %>%
  head(25) %>%
  rownames()

top_genes_mSEGs <- subsetted_mSEGs_as_features_loadings %>%
  arrange(desc(SumOfPCs)) %>%
  head(25) %>%
  rownames()

# Subset the Seurat objects based on the top 25 genes
subsetted_no_features_defined_subset <- subset(subsetted_no_features_defined, features = top_genes_no_features)
subsetted_mSEGs_as_features_subset <- subset(subsetted_mSEGs_as_features, features = top_genes_mSEGs)

# Extract the log-normalized counts
lognorm_counts_no_features <- as.matrix(GetAssayData(subsetted_no_features_defined_subset, slot = "data"))
lognorm_counts_mSEGs <- as.matrix(GetAssayData(subsetted_mSEGs_as_features_subset, slot = "data"))

# Normalize columns using the provided function
normalize_columns <- function(x) {
  (x - mean(x)) / sd(x)
}

lognorm_counts_no_features <- t(apply(lognorm_counts_no_features, 1, normalize_columns))
lognorm_counts_mSEGs <- t(apply(lognorm_counts_mSEGs, 1, normalize_columns))

# Transpose the matrices to have genes as columns and samples as rows
lognorm_counts_no_features <- t(lognorm_counts_no_features)
lognorm_counts_mSEGs <- t(lognorm_counts_mSEGs)

# Extract cell type annotations
cell_types_no_features <- as.factor(subsetted_no_features_defined_subset@meta.data$Type)
cell_types_mSEGs <- as.factor(subsetted_mSEGs_as_features_subset@meta.data$Type)

# Create a row annotation for cell types
row_annotation_no_features <- rowAnnotation(CellType = cell_types_no_features,
                                            col = list(CellType = c("2i" = "purple", "a2i" = "orange", "serum" = "darkgreen")))

row_annotation_mSEGs <- rowAnnotation(CellType = cell_types_mSEGs,
                                      col = list(CellType = c("2i" = "purple", "a2i" = "orange", "serum" = "darkgreen")))

# Create heatmaps
# Heatmap for subsetted_no_features_defined_subset
Heatmap(lognorm_counts_no_features, name = "LogNormalized Counts No Features",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_column_names = TRUE,
        heatmap_legend_param = list(title = "LogNormalized Counts"),
        column_title = "Genes", row_title = "Samples",
        right_annotation = row_annotation_no_features)

# Heatmap for subsetted_mSEGs_as_features_subset
Heatmap(lognorm_counts_mSEGs, name = "LogNormalized Counts mSEGs",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_column_names = TRUE,
        heatmap_legend_param = list(title = "LogNormalized Counts"),
        column_title = "Genes", row_title = "Samples",
        right_annotation = row_annotation_mSEGs)