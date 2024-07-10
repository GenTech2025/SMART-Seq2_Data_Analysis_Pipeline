# The filtered SCE object may contain annotation for the retained cells only
cell_type_fil <- colData(sce_raw_fil)$Type # Contains 515 cells 

# Similarly it can be extracted from seurat object
cell_type_fil_so <- so_raw_fil@meta.data$Type # Contains 515 cells 

# Extract the loadings
loading1 <- as.data.frame(Loadings(subsetted_mSEGs_as_features))
loading2 <- as.data.frame(Loadings(subsetted_no_features_defined))

# Extract the embeddings
embedding1 <- as.data.frame(Embeddings(subsetted_mSEGs_as_features))
embedding2 <- as.data.frame(Embeddings(subsetted_no_features_defined))

# Now we can calculate a third column which is the sum of absolute value of PC1 and PC2
loading1$combined <- abs(loading1$PC_1) + abs(loading1$PC_2)
loading2$combined <- abs(loading2$PC_1) + abs(loading2$PC_2)

# Extract the top genes based on the combined values
top_genes_1 <- loading1[order(loading1$combined, decreasing = TRUE), ][1:25,] # Top 25 genes where features have been explicitly defined
top_genes_2 <- loading2[order(loading2$combined, decreasing = TRUE), ]

# Just keep three columns 
top_genes_1_cln <- top_genes_1[,c("PC_1","PC_2","combined")] # Explicit features defined
top_genes_2_cln <- top_genes_2[,c("PC_1","PC_2","combined")] # No features defined

# Check the number of unique genes
unique_top_genes <- unique(c(rownames(top_genes_1_cln), rownames(top_genes_2_cln))) # There is a total of 48 unique genes 

# Map each gene to its origin PCA plot and extract combined values
gene_info <- lapply(unique_top_genes, function(gene) {
  origin <- c()
  combined_values <- c()
  if (gene %in% rownames(top_genes_1_cln)) {
    origin <- c(origin, "PCA Plot 1")
    combined_values <- c(combined_values, top_genes_1_cln[gene, "combined"])
  }
  if (gene %in% rownames(top_genes_2_cln)) {
    origin <- c(origin, "PCA Plot 2")
    combined_values <- c(combined_values, top_genes_2_cln[gene, "combined"])
  }
  combined_value <- mean(combined_values, na.rm = TRUE) # Average the combined values
  list(Gene = gene, Origin = paste(origin, collapse = ", "), Combined = combined_value)
})

# Convert the list to a data frame
top_genes <- do.call(rbind, lapply(gene_info, function(x) data.frame(Gene = x$Gene, Origin = x$Origin, Combined = x$Combined, stringsAsFactors = FALSE)))


################### DO THE SAME WITH EMBEDDINGS #########################

# Assuming embedding1 and embedding2 are data frames or matrices with embeddings as rows

# Calculate the combined values for embeddings
embedding1$combined <- abs(embedding1$PC_1) + abs(embedding1$PC_2)
embedding2$combined <- abs(embedding2$PC_1) + abs(embedding2$PC_2)

# Extract the top embeddings based on the combined values
top_embeddings_1 <- embedding1[order(embedding1$combined, decreasing = TRUE), ] # Top 25 embeddings for the first dataset
top_embeddings_2 <- embedding2[order(embedding2$combined, decreasing = TRUE), ]

# Just keep three columns 
top_embeddings_1_cln <- top_embeddings_1[,c("PC_1","PC_2","combined")] # First dataset
top_embeddings_2_cln <- top_embeddings_2[,c("PC_1","PC_2","combined")] # Second dataset

# Check the number of unique embeddings
unique_top_embeddings <- unique(c(rownames(top_embeddings_1_cln), rownames(top_embeddings_2_cln))) # Get unique embeddings

# Map each embedding to its origin PCA plot and extract combined values
embedding_info <- lapply(unique_top_embeddings, function(embedding) {
  origin <- c()
  combined_values <- c()
  if (embedding %in% rownames(top_embeddings_1_cln)) {
    origin <- c(origin, "PCA Plot 1")
    combined_values <- c(combined_values, top_embeddings_1_cln[embedding, "combined"])
  }
  if (embedding %in% rownames(top_embeddings_2_cln)) {
    origin <- c(origin, "PCA Plot 2")
    combined_values <- c(combined_values, top_embeddings_2_cln[embedding, "combined"])
  }
  combined_value <- mean(combined_values, na.rm = TRUE) # Average the combined values
  list(Embedding = embedding, Origin = paste(origin, collapse = ", "), Combined = combined_value)
})

# Convert the list to a data frame
top_embeddings <- do.call(rbind, lapply(embedding_info, function(x) data.frame(Embedding = x$Embedding, Origin = x$Origin, Combined = x$Combined, stringsAsFactors = FALSE)))

#### NOW I HAVE TWO DATAFRAMES - top_genes and top_embeddings

top_loadings <- top_genes

rownames(top_loadings) <- top_genes$Gene

top_loadings <- top_loadings$Combined
