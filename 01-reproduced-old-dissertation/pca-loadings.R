# Extract the loadings for all features used in PCA
pca_loading <- Loadings(so_decon_norm_seg, reduction = "pca")
# Convert to dataframe
pca_loading_df <- as.data.frame(pca_loading)
# Add a new column for ensembl ID
pca_loading_df$ensemblID <- rownames(pca_loading_df)

write.csv(pca_loading_df, file = "./pca_load_features.csv", row.names = FALSE)