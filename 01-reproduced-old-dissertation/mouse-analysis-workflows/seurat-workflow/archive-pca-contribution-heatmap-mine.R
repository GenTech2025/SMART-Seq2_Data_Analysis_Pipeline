### TEST GROUND

```{r fig.height=10, fig.width=12}
pca_embeddings <- Embeddings(subsetted_no_features_defined, reduction = "pca")

pca_loadings <- as.data.frame(subsetted_no_features_defined@reductions$pca@feature.loadings)

contribution_matrix <- pca_embeddings[, c("PC_1", "PC_2")] %*% t(pca_loadings[, c("PC_1", "PC_2")])

library(ComplexHeatmap)

# Create the heatmap
Heatmap(contribution_matrix, name = "Contribution",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_column_names = TRUE,
        heatmap_legend_param = list(title = "Contribution"))
```

#### USING Non-ORMALIZED DATA
```{r fig.height=15, fig.width=15}
pca_embeddings <- Embeddings(subsetted_no_features_defined, reduction = "pca")

pca_loadings <- as.data.frame(subsetted_no_features_defined@reductions$pca@feature.loadings)

contribution_matrix <- pca_embeddings[, c("PC_1", "PC_2")] %*% t(pca_loadings[, c("PC_1", "PC_2")])

cell_types <- subsetted_no_features_defined@meta.data$Type

# Function to normalize each column
normalize_columns <- function(x) {
  (x - mean(x)) / sd(x)
}

# Apply normalization to each column
normalized_matrix <- apply(contribution_matrix, 2, normalize_columns)

# Create a data frame for annotations
annotation_df <- data.frame(CellType = cell_types)
rownames(annotation_df) <- rownames(pca_embeddings)

# Create the heatmap with row annotations
Heatmap(normalized_matrix, name = "Contribution",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_column_names = TRUE,
        heatmap_legend_param = list(title = "Contribution"),
        right_annotation = rowAnnotation(CellType = annotation_df$CellType,
                                         col = list(CellType = c("a2i" = "orange", "2i" = "purple", "serum" = "green"))))
```
```{r fig.height=8, fig.width=10}
# Center the contribution matrix by columns
centered_contribution_matrix <- scale(contribution_matrix, center = TRUE, scale = FALSE)

Heatmap(contribution_matrix, name = "Contribution",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_column_names = TRUE,
        heatmap_legend_param = list(title = "Contribution"),
        right_annotation = rowAnnotation(CellType = annotation_df$CellType,
                                         col = list(CellType = c("a2i" = "orange", "2i" = "purple", "serum" = "green"))))
```


### USING NORMALIZED DATA
```{r}
# Extract PCA embeddings and loadings
pca_embeddings <- Embeddings(subsetted_no_features_defined, reduction = "pca")
pca_loadings <- as.data.frame(Loadings(subsetted_no_features_defined, reduction = "pca"))

# Compute contribution matrix
contribution_matrix <- pca_embeddings[, c("PC_1", "PC_2")] %*% t(pca_loadings[, c("PC_1", "PC_2")])

# Normalize contribution matrix
normalize_columns <- function(x) {
  (x - mean(x)) / sd(x)
}
normalized_matrix <- apply(contribution_matrix, 2, normalize_columns)

# Create annotation data frame
cell_types <- subsetted_no_features_defined@meta.data$Type
annotation_df <- data.frame(CellType = cell_types)
rownames(annotation_df) <- rownames(pca_embeddings)

# Generate heatmap
png("./heatmap_3.png", width = 12, height = 10, units = "in", res = 500)
Heatmap(normalized_matrix, name = "Normalized Contribution",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_column_names = TRUE,
        heatmap_legend_param = list(title = "Contribution"),
        right_annotation = rowAnnotation(CellType = annotation_df$CellType,
                                         col = list(CellType = c("a2i" = "orange", "2i" = "purple", "serum" = "green"))))
dev.off()
```


### TRYING TO DO SIMONS HEATMAP

```{r}
# Extract embeddings and loadings for both datasets
embeddings_1 <- Embeddings(subsetted_no_features_defined, reduction = "pca")
loadings_1 <- as.data.frame(subsetted_no_features_defined@reductions$pca@feature.loadings)

embeddings_2 <- Embeddings(subsetted_mSEGs_as_features, reduction = "pca")
loadings_2 <- as.data.frame(subsetted_mSEGs_as_features@reductions$pca@feature.loadings)

```


```{r}
# Get top 10 genes for PC1 and PC2 for both PCA plots
top_genes_1 <- unique(c(
  rownames(loadings_1[order(abs(loadings_1$PC_1), decreasing = TRUE),])[1:10],
  rownames(loadings_1[order(abs(loadings_1$PC_2), decreasing = TRUE),])[1:10]
))
top_genes_2 <- unique(c(
  rownames(loadings_2[order(abs(loadings_2$PC_1), decreasing = TRUE),])[1:10],
  rownames(loadings_2[order(abs(loadings_2$PC_2), decreasing = TRUE),])[1:10]
))

# Combine top genes and create a data frame with their origins
all_genes <- unique(c(top_genes_1, top_genes_2))

# Map each gene to its origin PCA plot
gene_origin <- sapply(all_genes, function(gene) {
  origin <- c()
  if (gene %in% top_genes_1) origin <- c(origin, "PCA Plot 1")
  if (gene %in% top_genes_2) origin <- c(origin, "PCA Plot 2")
  paste(origin, collapse = ", ")
})

top_genes <- data.frame(Gene = all_genes, Origin = gene_origin)

```

```{r}
# Calculate contribution matrices for both PCA plots using the combined list of top genes
valid_genes_1 <- intersect(top_genes$Gene, rownames(loadings_1))
valid_genes_2 <- intersect(top_genes$Gene, rownames(loadings_2))

contribution_matrix_1 <- embeddings_1[, c("PC_1", "PC_2")] %*% t(loadings_1[valid_genes_1, c("PC_1", "PC_2")])
contribution_matrix_2 <- embeddings_2[, c("PC_1", "PC_2")] %*% t(loadings_2[valid_genes_2, c("PC_1", "PC_2")])

```

```{r}
# Function to normalize each column
normalize_columns <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Apply normalization
normalized_matrix_1 <- apply(contribution_matrix_1, 2, normalize_columns)
normalized_matrix_2 <- apply(contribution_matrix_2, 2, normalize_columns)

```

```{r}
# Combine the normalized matrices into a single matrix
combined_matrix <- cbind(normalized_matrix_1[, valid_genes_1, drop = FALSE], normalized_matrix_2[, valid_genes_2, drop = FALSE])

# Ensure the row names (sample identifiers) are unique
rownames(combined_matrix) <- make.unique(rownames(combined_matrix))

# Create annotation data frame for cell types
cell_types_1 <- subsetted_no_features_defined@meta.data$Type
cell_types_2 <- subsetted_mSEGs_as_features@meta.data$Type

cell_types <- c(cell_types_1, cell_types_2)

# Ensure the annotation data frame matches the combined matrix row names
annotation_df <- data.frame(CellType = cell_types)
rownames(annotation_df) <- rownames(combined_matrix)

# Generate heatmap
library(ComplexHeatmap)

png("./combined_heatmap.png", width = 12, height = 10, units = "in", res = 500)
Heatmap(combined_matrix, name = "Normalized Contribution",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = TRUE, show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "Contribution"),
        top_annotation = HeatmapAnnotation(
          PlotOrigin = anno_simple(top_genes$Origin, col = c("PCA Plot 1" = "blue", "PCA Plot 2" = "red", "PCA Plot 1, PCA Plot 2" = "purple"))
        ),
        right_annotation = rowAnnotation(CellType = annotation_df$CellType,
                                         col = list(CellType = c("a2i" = "orange", "2i" = "purple", "serum" = "green"))))
dev.off()


```



### REVISED CODE
```{r}
# Get top 10 genes for PC1 and PC2 for both PCA plots
top_genes_1 <- unique(c(
  rownames(loadings_1[order(abs(loadings_1$PC_1), decreasing = TRUE),])[1:10],
  rownames(loadings_1[order(abs(loadings_1$PC_2), decreasing = TRUE),])[1:10]
))
top_genes_2 <- unique(c(
  rownames(loadings_2[order(abs(loadings_2$PC_1), decreasing = TRUE),])[1:10],
  rownames(loadings_2[order(abs(loadings_2$PC_2), decreasing = TRUE),])[1:10]
))

# Combine top genes and create a data frame with their origins
all_genes <- unique(c(top_genes_1, top_genes_2))

# Map each gene to its origin PCA plot
gene_origin <- sapply(all_genes, function(gene) {
  origin <- c()
  if (gene %in% top_genes_1) origin <- c(origin, "PCA Plot 1")
  if (gene %in% top_genes_2) origin <- c(origin, "PCA Plot 2")
  paste(origin, collapse = ", ")
})

top_genes <- data.frame(Gene = all_genes, Origin = gene_origin, stringsAsFactors = FALSE)

# Align matrices to have the same set of genes
valid_genes <- intersect(top_genes$Gene, intersect(rownames(loadings_1), rownames(loadings_2)))

# Recreate top_genes to ensure only valid genes are included
top_genes <- top_genes[top_genes$Gene %in% valid_genes, ]

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

# Ensure both normalized matrices have the same columns
colnames(normalized_matrix_1) <- valid_genes
colnames(normalized_matrix_2) <- valid_genes

# Combine the normalized matrices into a single matrix
combined_matrix <- rbind(normalized_matrix_1, normalized_matrix_2)

# Ensure the row names (sample identifiers) are unique
rownames(combined_matrix) <- make.unique(rownames(combined_matrix))

# Verify lengths
print(paste("Number of rows in combined_matrix:", nrow(combined_matrix)))
print(paste("Number of cell types:", length(c(subsetted_no_features_defined@meta.data$Type, subsetted_mSEGs_as_features@meta.data$Type))))
print(paste("Number of columns in combined_matrix:", ncol(combined_matrix)))
print(paste("Number of top genes:", nrow(top_genes)))

# Ensure top_genes are aligned with columns
if (ncol(combined_matrix) != nrow(top_genes)) {
  stop("Mismatch between number of columns in combined_matrix and number of top genes")
}

# Create annotation data frame for cell types
cell_types_1 <- subsetted_no_features_defined@meta.data$Type
cell_types_2 <- subsetted_mSEGs_as_features@meta.data$Type

cell_types <- c(cell_types_1, cell_types_2)

# Ensure the lengths match
if (nrow(combined_matrix) != length(cell_types)) {
  stop("Mismatch between number of rows in combined_matrix and length of cell_types")
}

annotation_df <- data.frame(CellType = cell_types, stringsAsFactors = FALSE)
rownames(annotation_df) <- rownames(combined_matrix)

# Generate heatmap
library(ComplexHeatmap)

# Ensure top_genes are ordered correctly
top_genes <- top_genes[match(colnames(combined_matrix), top_genes$Gene), ]

png("./combined_heatmap.png", width = 12, height = 10, units = "in", res = 500)
Heatmap(combined_matrix, name = "Normalized Contribution",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "Contribution"),
        top_annotation = HeatmapAnnotation(
          PlotOrigin = anno_simple(top_genes$Origin, col = c("PCA Plot 1" = "blue", "PCA Plot 2" = "red", "PCA Plot 1, PCA Plot 2" = "purple"))
        ),
        right_annotation = rowAnnotation(CellType = annotation_df$CellType,
                                         col = list(CellType = c("a2i" = "orange", "2i" = "purple", "serum" = "green"))))
dev.off()

```



### Corrected code
```{r}
# Extract embeddings and loadings for both datasets
embeddings_1 <- Embeddings(subsetted_no_features_defined, reduction = "pca")
loadings_1 <- as.data.frame(subsetted_no_features_defined@reductions$pca@feature.loadings)

embeddings_2 <- Embeddings(subsetted_mSEGs_as_features, reduction = "pca")
loadings_2 <- as.data.frame(subsetted_mSEGs_as_features@reductions$pca@feature.loadings)

########################################################################################

# Calculate combined contribution for PC1 and PC2
loadings_1$combined <- abs(loadings_1$PC_1) + abs(loadings_1$PC_2)
loadings_2$combined <- abs(loadings_2$PC_1) + abs(loadings_2$PC_2)

# Get top 10 genes based on the combined contribution
top_genes_1 <- rownames(loadings_1[order(loadings_1$combined, decreasing = TRUE), ])[1:10]
top_genes_2 <- rownames(loadings_2[order(loadings_2$combined, decreasing = TRUE), ])[1:10]

# Combine top genes and create a data frame with their origins
all_genes <- unique(c(top_genes_1, top_genes_2))

# Map each gene to its origin PCA plot
gene_origin <- sapply(all_genes, function(gene) {
  origin <- c()
  if (gene %in% top_genes_1) origin <- c(origin, "PCA Plot 1")
  if (gene %in% top_genes_2) origin <- c(origin, "PCA Plot 2")
  paste(origin, collapse = ", ")
})

top_genes <- data.frame(Gene = all_genes, Origin = gene_origin, stringsAsFactors = FALSE)


# Align matrices to have the same set of genes
valid_genes <- intersect(top_genes$Gene, intersect(rownames(loadings_1), rownames(loadings_2)))

# Recreate top_genes to ensure only valid genes are included
top_genes <- top_genes[top_genes$Gene %in% valid_genes, ]

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


# Ensure both normalized matrices have the same columns
colnames(normalized_matrix_1) <- valid_genes
colnames(normalized_matrix_2) <- valid_genes

# Combine the normalized matrices into a single matrix
combined_matrix <- rbind(normalized_matrix_1, normalized_matrix_2)

# Ensure the row names (sample identifiers) are unique
rownames(combined_matrix) <- make.unique(c(rownames(normalized_matrix_1), rownames(normalized_matrix_2)))

# Verify lengths
print(paste("Number of rows in combined_matrix:", nrow(combined_matrix)))
print(paste("Number of cell types:", length(c(subsetted_no_features_defined@meta.data$Type, subsetted_mSEGs_as_features@meta.data$Type))))
print(paste("Number of columns in combined_matrix:", ncol(combined_matrix)))
print(paste("Number of top genes:", nrow(top_genes)))

# Ensure top_genes are aligned with columns
if (ncol(combined_matrix) != nrow(top_genes)) {
  stop("Mismatch between number of columns in combined_matrix and number of top genes")
}

# Create annotation data frame for cell types
cell_types_1 <- subsetted_no_features_defined@meta.data$Type
cell_types_2 <- subsetted_mSEGs_as_features@meta.data$Type

cell_types <- c(cell_types_1, cell_types_2)

# Ensure the lengths match
if (nrow(combined_matrix) != length(cell_types)) {
  stop("Mismatch between number of rows in combined_matrix and length of cell_types")
}

annotation_df <- data.frame(CellType = cell_types, stringsAsFactors = FALSE)
rownames(annotation_df) <- rownames(combined_matrix)

# Generate heatmap
library(ComplexHeatmap)

# Ensure top_genes are ordered correctly
top_genes <- top_genes[match(colnames(combined_matrix), top_genes$Gene), ]

Heatmap(combined_matrix, name = "Normalized Contribution",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "Contribution"),
        top_annotation = HeatmapAnnotation(
          PlotOrigin = anno_simple(top_genes$Origin, col = c("PCA Plot 1" = "blue", "PCA Plot 2" = "red", "PCA Plot 1, PCA Plot 2" = "purple"))
        ),
        right_annotation = rowAnnotation(CellType = annotation_df$CellType,
                                         col = list(CellType = c("a2i" = "orange", "2i" = "purple", "serum" = "green"))))
```





---
  
  ####JKJKJKJKJKJKJKJKJKJKJKKKJKJKJKJKKKJKJKJKKJKJKK
  
  # Try to make seperate heatmaps for both PCA plots
  #### Extract the PCA loadings and embeddings
  ```{r}
# Subsetted Seurat object with no explicit feature definition
subsetted_no_features_defined_loadings <- Loadings(subsetted_no_features_defined, reduction = "pca")[,1:2]
subsetted_no_features_defined_embeddings <- Embeddings(subsetted_no_features_defined, reduction = "pca")[,1:2]

# Subsetted Seurat object with explicit definition of feature
subsetted_mSEGs_as_features_loadings <- Loadings(subsetted_mSEGs_as_features, reduction = "pca")[,1:2]
subsetted_mSEGs_as_features_embeddings <- Embeddings(subsetted_mSEGs_as_features, reduction = "pca")[,1:2]
```

#### Convert to dataframes
```{r}
# Subsetted Seurat object with no explicit feature definition
subsetted_no_features_defined_loadings <- as.data.frame(subsetted_no_features_defined_loadings)
subsetted_no_features_defined_embeddings <- as.data.frame(subsetted_no_features_defined_embeddings)
# Subsetted Seurat object with explicit definition of feature
subsetted_mSEGs_as_features_loadings <- as.data.frame(subsetted_mSEGs_as_features_loadings)
subsetted_mSEGs_as_features_embeddings <- as.data.frame(subsetted_mSEGs_as_features_embeddings)
```

#### Add a new column to the loading dataframes which is the sum of first two principal components
```{r}
# Subsetted Seurat object with no explicit feature definition
subsetted_no_features_defined_loadings$SumOfPCs <- abs(subsetted_no_features_defined_loadings$PC_1) + abs(subsetted_no_features_defined_loadings$PC_2)
# Subsetted Seurat object with explicit definition of feature
subsetted_mSEGs_as_features_loadings$SumofPCs <- abs(subsetted_mSEGs_as_features_loadings)$PC_1 + abs(subsetted_mSEGs_as_features_loadings)$PC_2
```

#### Add a new column to the embedding dataframes which is the sum of first two principal components
```{r}
subsetted_no_features_defined_embeddings$SumOfPCs <- abs(subsetted_no_features_defined_embeddings$PC_1) + abs(subsetted_no_features_defined_embeddings$PC_2)

subsetted_mSEGs_as_features_embeddings$SumOfPCs <- abs(subsetted_mSEGs_as_features_embeddings$PC_1) + abs(subsetted_mSEGs_as_features_embeddings$PC_2)
```

#### Arrange the dataframes in descending order of the new column values
```{r}
# Arrange all the dataframes in descending order of SumOfPCs values
subsetted_no_features_defined_loadings <- subsetted_no_features_defined_loadings %>%
  arrange(desc(SumOfPCs))

subsetted_mSEGs_as_features_loadings <- subsetted_mSEGs_as_features_loadings %>%
  arrange(desc(SumOfPCs))

subsetted_no_features_defined_embeddings <- subsetted_no_features_defined_embeddings %>%
  arrange(desc(SumOfPCs))

subsetted_mSEGs_as_features_embeddings <- subsetted_mSEGs_as_features_embeddings %>%
  arrange(desc(SumOfPCs))
```

## Compute the contribution matrix of Seurat Normalized + Subsetted + No Features Defined
```{r}
# Ensure the relevant columns are numeric
subsetted_no_features_defined_embeddings <- subsetted_no_features_defined_embeddings %>%
  mutate(PC_1 = as.numeric(PC_1), PC_2 = as.numeric(PC_2))

subsetted_no_features_defined_loadings <- subsetted_no_features_defined_loadings %>%
  mutate(PC_1 = as.numeric(PC_1), PC_2 = as.numeric(PC_2))

# Convert the selected columns to matrices
embeddings_matrix <- as.matrix(subsetted_no_features_defined_embeddings[, c("PC_1", "PC_2")])
loadings_matrix <- as.matrix(subsetted_no_features_defined_loadings[, c("PC_1", "PC_2")])

# Perform matrix multiplication
contribution_matrix_1 <- embeddings_matrix %*% t(loadings_matrix)

# Function to normalize each column
normalize_columns <- function(x) {
  (x - mean(x)) / sd(x)
}

# Apply normalization to each column
normalized_matrix <- apply(contribution_matrix_1, 2, normalize_columns)
```

## Create Heatmap to visualize the matrix
```{r fig.height=10, fig.width=12}
# Extract the cell type annotation
cell_types <- subsetted_no_features_defined@meta.data$Type

# Create a data frame for annotations
annotation_df <- data.frame(CellType = cell_types)
rownames(annotation_df) <- rownames(subsetted_no_features_defined_embeddings)

# Create the heatmap with row annotations
Heatmap(normalized_matrix, name = "Contribution",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_column_names = TRUE,
        heatmap_legend_param = list(title = "Contribution"),
        right_annotation = rowAnnotation(CellType = annotation_df$CellType,
                                         col = list(CellType = c("a2i" = "orange", "2i" = "purple", "serum" = "green"))))
```

## DO the same for the other PCA Plot
```{r fig.height=10, fig.width=12}
## Compute the contribution matrix of Seurat Normalized + Subsetted + mSEGs as Features
# Ensure the relevant columns are numeric
subsetted_mSEGs_as_features_embeddings <- subsetted_mSEGs_as_features_embeddings %>%
  mutate(PC_1 = as.numeric(PC_1), PC_2 = as.numeric(PC_2))

subsetted_mSEGs_as_features_loadings <- subsetted_mSEGs_as_features_loadings %>%
  mutate(PC_1 = as.numeric(PC_1), PC_2 = as.numeric(PC_2))

# Convert the selected columns to matrices
embeddings_matrix <- as.matrix(subsetted_mSEGs_as_features_embeddings[, c("PC_1", "PC_2")])
loadings_matrix <- as.matrix(subsetted_mSEGs_as_features_loadings[, c("PC_1", "PC_2")])

# Use only the top 25 rows of the loading matrix
loadings_matrix_top25 <- loadings_matrix[1:25, ]

# Perform matrix multiplication
contribution_matrix_2 <- embeddings_matrix %*% t(loadings_matrix_top25)

# Function to normalize each column
normalize_columns <- function(x) {
  (x - mean(x)) / sd(x)
}

# Apply normalization to each column
normalized_matrix <- apply(contribution_matrix_2, 2, normalize_columns)

## Create Heatmap to visualize the matrix
# Extract the cell type annotation
cell_types <- subsetted_mSEGs_as_features@meta.data$Type

# Create a data frame for annotations
annotation_df <- data.frame(CellType = cell_types)
rownames(annotation_df) <- rownames(subsetted_mSEGs_as_features_embeddings)

# Create the heatmap with row annotations
Heatmap(normalized_matrix, name = "Contribution",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_column_names = TRUE,
        heatmap_legend_param = list(title = "Contribution"),
        right_annotation = rowAnnotation(CellType = annotation_df$CellType,
                                         col = list(CellType = c("a2i" = "orange", "2i" = "purple", "serum" = "green"))))
```

# Try to combine the contribution data frames
```{r fig.height=12, fig.height=14}
# Function to normalize each column
normalize_columns <- function(x) {
  (x - mean(x)) / sd(x)
}

# Ensure the relevant columns are numeric for the first Seurat object
subsetted_no_features_defined_embeddings <- subsetted_no_features_defined_embeddings %>%
  mutate(PC_1 = as.numeric(PC_1), PC_2 = as.numeric(PC_2))

subsetted_no_features_defined_loadings <- subsetted_no_features_defined_loadings %>%
  mutate(PC_1 = as.numeric(PC_1), PC_2 = as.numeric(PC_2))

# Convert the selected columns to matrices for the first Seurat object
embeddings_matrix_1 <- as.matrix(subsetted_no_features_defined_embeddings[, c("PC_1", "PC_2")])
loadings_matrix_1 <- as.matrix(subsetted_no_features_defined_loadings[, c("PC_1", "PC_2")])

# Perform matrix multiplication for the first Seurat object
contribution_matrix_1 <- embeddings_matrix_1 %*% t(loadings_matrix_1)

# Normalize the contribution matrix for the first Seurat object
normalized_matrix_1 <- apply(contribution_matrix_1, 2, normalize_columns)

# Ensure the relevant columns are numeric for the second Seurat object
subsetted_mSEGs_as_features_embeddings <- subsetted_mSEGs_as_features_embeddings %>%
  mutate(PC_1 = as.numeric(PC_1), PC_2 = as.numeric(PC_2))

subsetted_mSEGs_as_features_loadings <- subsetted_mSEGs_as_features_loadings %>%
  mutate(PC_1 = as.numeric(PC_1), PC_2 = as.numeric(PC_2))

# Convert the selected columns to matrices for the second Seurat object
embeddings_matrix_2 <- as.matrix(subsetted_mSEGs_as_features_embeddings[, c("PC_1", "PC_2")])
loadings_matrix_2 <- as.matrix(subsetted_mSEGs_as_features_loadings[, c("PC_1", "PC_2")])

# Use only the top 25 rows of the loading matrix for the second Seurat object
loadings_matrix_top25 <- loadings_matrix_2[1:25, ]

# Perform matrix multiplication for the second Seurat object
contribution_matrix_2 <- embeddings_matrix_2 %*% t(loadings_matrix_top25)

# Normalize the contribution matrix for the second Seurat object
normalized_matrix_2 <- apply(contribution_matrix_2, 2, normalize_columns)

# Convert both normalized matrices to data frames and add suffixes to the column names
df1 <- as.data.frame(normalized_matrix_1)
colnames(df1) <- paste0(colnames(df1), "_obj1")

df2 <- as.data.frame(normalized_matrix_2)
colnames(df2) <- paste0(colnames(df2), "_obj2")

# Combine the data frames
combined_df <- cbind(df1, df2)

# Handle common genes (columns) in both data frames
common_genes <- intersect(colnames(df1), colnames(df2))
for (gene in common_genes) {
  combined_df[[gene]] <- df1[[gene]] + df2[[gene]]
}

# Create the heatmap with the combined data frame
Heatmap(combined_df, name = "Combined Contribution",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_column_names = TRUE,
        heatmap_legend_param = list(title = "Combined Contribution"),
        right_annotation = rowAnnotation(CellType = annotation_df$CellType,
                                         col = list(CellType = c("a2i" = "orange", "2i" = "purple", "serum" = "green"))))
```

```{r fig.height=15, fig.width=15}
# Create a vector for gene origins
gene_origins <- c(rep("obj1", ncol(df1)), rep("obj2", ncol(df2)))
gene_origins[which(names(combined_df) %in% common_genes)] <- "both"

# Create a column annotation for gene origins
column_annotation <- HeatmapAnnotation(
  GeneOrigin = gene_origins,
  col = list(GeneOrigin = c("obj1" = "blue", "obj2" = "red", "both" = "purple"))
)

# Create the heatmap with the combined data frame and annotations
Heatmap(combined_df, name = "Combined Contribution",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_column_names = TRUE,
        heatmap_legend_param = list(title = "Combined Contribution"),
        top_annotation = column_annotation,
        right_annotation = rowAnnotation(CellType = annotation_df$CellType,
                                         col = list(CellType = c("a2i" = "orange", "2i" = "purple", "serum" = "green"))))
```


#######################################################################################################################
library(Seurat)
library(scMerge)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(corrplot)
library(scater)
library(scran)
library(ComplexHeatmap)
library(biomaRt)
library(dplyr)
library(knitr)

data_directory <- "/home/s2599932/study-datasets/mouse-dataset/E-MTAB-2600/single_es_tximport.rds"
annotation_directory <- "/home/s2599932/study-datasets/mouse-dataset/E-MTAB-2600/E-MTAB-2600.targets.txt"
output_files_directory <- "./"
output_plots_directory <- "./"

annotation <- read.table(annotation_directory, sep = "\t", header = TRUE, row.names = "ERR")
raw_data <- readRDS(data_directory)

data("segList_ensemblGeneID", package = "scMerge")
mSEG <- segList_ensemblGeneID$mouse$mouse_scSEG # Extracting SEGs identified in mouse by scRNA Seq

sce_raw <- SingleCellExperiment(assays=list(counts = raw_data$counts, abundance = raw_data$abundance))
colData(sce_raw) <- DataFrame(annotation)

is_spike <- grepl("^ERCC", rownames(sce_raw))
mitogenes <- "ENSMUSG00000064336|ENSMUSG00000064337|ENSMUSG00000064338|ENSMUSG00000064339|ENSMUSG00000064340|ENSMUSG00000064341|ENSMUSG00000064342|ENSMUSG00000064343|ENSMUSG00000064344|ENSMUSG00000064345|ENSMUSG00000064346|ENSMUSG00000064347|ENSMUSG00000064348|ENSMUSG00000064349|ENSMUSG00000064350|ENSMUSG00000064351|ENSMUSG00000064352|ENSMUSG00000064353|ENSMUSG00000064354|ENSMUSG00000064355|ENSMUSG00000064356|ENSMUSG00000064357|ENSMUSG00000064358|ENSMUSG00000064359|ENSMUSG00000064360|ENSMUSG00000064361|ENSMUSG00000065947|ENSMUSG00000064363|ENSMUSG00000064364|ENSMUSG00000064365|ENSMUSG00000064366|ENSMUSG00000064367|ENSMUSG00000064368|ENSMUSG00000064369|ENSMUSG00000064370|ENSMUSG00000064371|ENSMUSG00000064372"
is_mito <- grepl(mitogenes, rownames(sce_raw))

sce_raw$qc_stats <- perCellQCMetrics(sce_raw, subsets=list(ERCC=is_spike, Mt=is_mito))

# Calculating the filtering metrics
libsize_drop <- isOutlier(sce_raw$qc_stats$sum, nmads = 3, type = "lower", log = TRUE)
feature_drop <- isOutlier(sce_raw$qc_stats$detected, nmads = 3, type = "lower", log = TRUE)
mito_drop <- isOutlier(sce_raw$qc_stats$subsets_Mt_percent, nmads = 3, type = "higher")
spike_drop <- isOutlier(sce_raw$qc_stats$subsets_ERCC_percent, nmads = 3, type = "higher")

sce_raw_copy <- sce_raw # Save a non-filtered copy

# Filtering the SingleCellExperiment (SCE) object
sce_raw_fil <- sce_raw[, !(libsize_drop | feature_drop | mito_drop | spike_drop)]
sce_raw_fil_copy <- sce_raw_fil # Save a filtered copy

# Save the filtered cells
filtered_cells <- colnames(sce_raw_fil)

# Create Seurat Object
raw_counts <- as.data.frame(raw_data$counts)

so_raw <- CreateSeuratObject(counts = raw_counts)
so_raw <- AddMetaData(so_raw, metadata = annotation)

so_raw_copy <- so_raw

# Only include samples that passed QC
so_raw_fil <- so_raw[,filtered_cells]
so_raw_fil_copy <- so_raw_fil

# Normalize and Scale the data
so_norm_fil <- NormalizeData(so_raw_fil, verbose = FALSE)

so_norm_fil <- FindVariableFeatures(so_norm_fil, selection.method = "vst", verbose = FALSE)

so_norm_fil <- ScaleData(so_norm_fil,features = rownames(so_norm_fil), verbose = FALSE)

# Select the cell types we are interested in
conditions_of_interest <- c("2i", "a2i", "serum")
cells_of_interest <- WhichCells(so_norm_fil, expression = Type %in% conditions_of_interest)

# Subset the seurat object so that it contains only those cell types
so_norm_fil <- subset(so_norm_fil, cells = cells_of_interest)

# Copy the seurat object for different combination of PCAs
non_subsetted_no_features_defined <- so_norm_fil
non_subsetted_all_genes_as_features <- so_norm_fil
non_subsetted_mSEGs_as_features <- so_norm_fil
subsetted_no_features_defined <- so_norm_fil
subsetted_mSEGs_as_features <- so_norm_fil

# Set Custom color for PCA plot
custom_colors <- c(
  "2i" = "#1f77b4",  # Blue
  "a2i" = "#2ca02c",  # Green
  "serum" = "#9467bd",   # Purple
)

# Select the outlying cells that needs to be removed
outliers_df <- read.csv("../sample_outliers_SN_NS_NF.csv", header = FALSE)
outliers <- outliers_df$X

# Generate the PCA plots with different combination of parameters (only two included in this code)

# subsetted no features defined
set.seed(194)

subsetted_no_features_defined <- subset(subsetted_no_features_defined, cells = setdiff(Cells(subsetted_no_features_defined), outliers))

subsetted_no_features_defined <- subset(subsetted_no_features_defined, features = mSEG)

subsetted_no_features_defined <- RunPCA(object = subsetted_no_features_defined, verbose = FALSE)

p4 <- DimPlot(subsetted_no_features_defined, reduction = "pca", group.by = "Type") +
  ggtitle("Seurat Normalized + Subsetted + No features defined (E-MTAB-2600) + No 2i outliers") +
  theme_minimal(base_size = 15) +  # Increase base font size for better readability
  scale_color_manual(values = custom_colors) +  # Use a specific viridis color option
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),  # Center and bold the title
    axis.title.x = element_text(size = 15),  # Increase x-axis title size
    axis.title.y = element_text(size = 15),  # Increase y-axis title size
    axis.text = element_text(size = 12),  # Increase axis text size
    legend.title = element_text(size = 12),  # Increase legend title size
    legend.text = element_text(size = 12)  # Increase legend text size
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # Increase legend point size


print(p4)

# subsetted mSEGs as features
set.seed(195)

subsetted_mSEGs_as_features <- subset(subsetted_mSEGs_as_features, cells = setdiff(Cells(subsetted_mSEGs_as_features), outliers))

subsetted_mSEGs_as_features <- subset(subsetted_mSEGs_as_features, features = mSEG)

subsetted_mSEGs_as_features <- RunPCA(object = subsetted_mSEGs_as_features, features = mSEG, verbose = FALSE)

p5 <- DimPlot(subsetted_mSEGs_as_features, reduction = "pca", group.by = "Type") +
  ggtitle("Seurat Normalized + Subsetted + mSEGs as features (E-MTAB-2600) + No 2i outliers") +
  theme_minimal(base_size = 15) +  # Increase base font size for better readability
  scale_color_manual(values = custom_colors) +  # Use a specific viridis color option
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),  # Center and bold the title
    axis.title.x = element_text(size = 15),  # Increase x-axis title size
    axis.title.y = element_text(size = 15),  # Increase y-axis title size
    axis.text = element_text(size = 12),  # Increase axis text size
    legend.title = element_text(size = 12),  # Increase legend title size
    legend.text = element_text(size = 12)  # Increase legend text size
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))  # Increase legend point size


print(p5)

# Now calculate the contribution matrix for the two PCA plots and then combine them

# Subsetted Seurat object with no explicit feature definition
subsetted_no_features_defined_loadings <- Loadings(subsetted_no_features_defined, reduction = "pca")[,1:2]
subsetted_no_features_defined_embeddings <- Embeddings(subsetted_no_features_defined, reduction = "pca")[,1:2]

# Subsetted Seurat object with explicit definition of feature
subsetted_mSEGs_as_features_loadings <- Loadings(subsetted_mSEGs_as_features, reduction = "pca")[,1:2]
subsetted_mSEGs_as_features_embeddings <- Embeddings(subsetted_mSEGs_as_features, reduction = "pca")[,1:2]

# Subsetted Seurat object with no explicit feature definition
subsetted_no_features_defined_loadings <- as.data.frame(subsetted_no_features_defined_loadings)
subsetted_no_features_defined_embeddings <- as.data.frame(subsetted_no_features_defined_embeddings)
# Subsetted Seurat object with explicit definition of feature
subsetted_mSEGs_as_features_loadings <- as.data.frame(subsetted_mSEGs_as_features_loadings)
subsetted_mSEGs_as_features_embeddings <- as.data.frame(subsetted_mSEGs_as_features_embeddings)

# Add the first two principal components and store them in a new column

# Subsetted Seurat object with no explicit feature definition
subsetted_no_features_defined_loadings$SumOfPCs <- abs(subsetted_no_features_defined_loadings$PC_1) + abs(subsetted_no_features_defined_loadings$PC_2)
subsetted_no_features_defined_embeddings$SumOfPCs <- abs(subsetted_no_features_defined_embeddings$PC_1) + abs(subsetted_no_features_defined_embeddings$PC_2)

# Subsetted Seurat object with explicit definition of feature
subsetted_mSEGs_as_features_loadings$SumofPCs <- abs(subsetted_mSEGs_as_features_loadings)$PC_1 + abs(subsetted_mSEGs_as_features_loadings)$PC_2
subsetted_mSEGs_as_features_embeddings$SumOfPCs <- abs(subsetted_mSEGs_as_features_embeddings$PC_1) + abs(subsetted_mSEGs_as_features_embeddings$PC_2)

# Arrange all the dataframes in descending order of SumOfPCs values
subsetted_no_features_defined_loadings <- subsetted_no_features_defined_loadings %>%
  arrange(desc(SumOfPCs))

subsetted_mSEGs_as_features_loadings <- subsetted_mSEGs_as_features_loadings %>%
  arrange(desc(SumOfPCs))

subsetted_no_features_defined_embeddings <- subsetted_no_features_defined_embeddings %>%
  arrange(desc(SumOfPCs))

subsetted_mSEGs_as_features_embeddings <- subsetted_mSEGs_as_features_embeddings %>%
  arrange(desc(SumOfPCs))

# Function to normalize each column
normalize_columns <- function(x) {
  (x - mean(x)) / sd(x)
}

# Ensure the relevant columns are numeric for the first Seurat object
subsetted_no_features_defined_embeddings <- subsetted_no_features_defined_embeddings %>%
  mutate(PC_1 = as.numeric(PC_1), PC_2 = as.numeric(PC_2))

subsetted_no_features_defined_loadings <- subsetted_no_features_defined_loadings %>%
  mutate(PC_1 = as.numeric(PC_1), PC_2 = as.numeric(PC_2))

# Convert the selected columns to matrices for the first Seurat object
embeddings_matrix_1 <- as.matrix(subsetted_no_features_defined_embeddings[, c("PC_1", "PC_2")])
loadings_matrix_1 <- as.matrix(subsetted_no_features_defined_loadings[, c("PC_1", "PC_2")])

# Perform matrix multiplication for the first Seurat object
contribution_matrix_1 <- embeddings_matrix_1 %*% t(loadings_matrix_1)

# Normalize the contribution matrix for the first Seurat object
normalized_matrix_1 <- apply(contribution_matrix_1, 2, normalize_columns)

# Ensure the relevant columns are numeric for the second Seurat object
subsetted_mSEGs_as_features_embeddings <- subsetted_mSEGs_as_features_embeddings %>%
  mutate(PC_1 = as.numeric(PC_1), PC_2 = as.numeric(PC_2))

subsetted_mSEGs_as_features_loadings <- subsetted_mSEGs_as_features_loadings %>%
  mutate(PC_1 = as.numeric(PC_1), PC_2 = as.numeric(PC_2))

# Convert the selected columns to matrices for the second Seurat object
embeddings_matrix_2 <- as.matrix(subsetted_mSEGs_as_features_embeddings[, c("PC_1", "PC_2")])
loadings_matrix_2 <- as.matrix(subsetted_mSEGs_as_features_loadings[, c("PC_1", "PC_2")])

# Use only the top 25 rows of the loading matrix for the second Seurat object
loadings_matrix_top25 <- loadings_matrix_2[1:25, ]

# Perform matrix multiplication for the second Seurat object
contribution_matrix_2 <- embeddings_matrix_2 %*% t(loadings_matrix_top25)

# Normalize the contribution matrix for the second Seurat object
normalized_matrix_2 <- apply(contribution_matrix_2, 2, normalize_columns)

# Convert both normalized matrices to data frames and add suffixes to the column names
df1 <- as.data.frame(normalized_matrix_1)
colnames(df1) <- paste0(colnames(df1), "_obj1")

df2 <- as.data.frame(normalized_matrix_2)
colnames(df2) <- paste0(colnames(df2), "_obj2")

# Combine the data frames
combined_df <- cbind(df1, df2)

# Handle common genes (columns) in both data frames
common_genes <- intersect(colnames(df1), colnames(df2))
for (gene in common_genes) {
  combined_df[[gene]] <- df1[[gene]] + df2[[gene]]
}

# Create the heatmap with the combined data frame
Heatmap(combined_df, name = "Combined Contribution",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = FALSE, show_column_names = TRUE,
        heatmap_legend_param = list(title = "Combined Contribution"),
        right_annotation = rowAnnotation(CellType = annotation_df$CellType,
                                         col = list(CellType = c("a2i" = "orange", "2i" = "purple", "serum" = "green"))))
#######################################################################################################################