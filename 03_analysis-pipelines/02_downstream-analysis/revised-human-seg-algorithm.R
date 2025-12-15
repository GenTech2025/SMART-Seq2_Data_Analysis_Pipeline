# Title: finalized-human-ranked-SEG-enrichment
# Author: Sourav Roy
# Date: 2024-08-01

# Load required libraries
library(tidyverse)
library(scMerge)
library(SingleCellExperiment)
library(scater)
library(scran)
library(biomaRt)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Seurat)
library(fgsea)
library(msigdb)
library(GSEABase)
library(jsonlite)
library(ggplot2)
library(pheatmap)

# Define File Paths
dataset_one_dir <- "/home/s2599932/study-datasets/human-dataset/E-MTAB-3929/E-MTAB-3929_counts.csv"
dataset_two_dir <- "/home/s2599932/study-datasets/human-dataset/GSE75748/GSE75748_counts.csv"
dataset_one_annotation_dir <- "/home/s2599932/study-datasets/human-dataset/E-MTAB-3929/E-MTAB-3929_targets.csv"
dataset_two_annotation_dir <- "/home/s2599932/study-datasets/human-dataset/GSE75748/GSE75748_targets.csv"

# Load the expression data and metadata
# Load data for the first set (E-MTAB-3929)
raw_one <- read.csv(dataset_one_dir, header = TRUE, row.names = 1)
annotation_one <- read.csv(dataset_one_annotation_dir, header = TRUE, row.names = "run_accession")
annotation_one <- annotation_one[match(colnames(raw_one), rownames(annotation_one)), ]

# Load data for the second set (GSE75748)
raw_two <- read.csv(dataset_two_dir, header = TRUE, row.names = 1)
annotation_two <- read.csv(dataset_two_annotation_dir, header = TRUE, row.names = "run_accession")
annotation_two <- annotation_two[match(colnames(raw_two), rownames(annotation_two)), ]

# Creating the single cell experiment object
sce_raw_one <- SingleCellExperiment(
  assays = list(counts = as.matrix(raw_one)),
  colData = annotation_one
)

sce_raw_two <- SingleCellExperiment(
  assays = list(counts = as.matrix(raw_two)),
  colData = annotation_two
)

# Extract the human stably expressed genes from scMerge package
data("segList_ensemblGeneID", package = "scMerge")
hSEG <- segList_ensemblGeneID$human$human_scSEG

# Quality Control
# Retrieve the mitochondrial genes for humans
# Connect to the Ensembl BioMart database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query for mitochondrial genes
mitogenes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "chromosome_name",
  values = "MT",
  mart = ensembl
)

# Extract the ensembl_gene_id column and combine into a single string
mitogenes_ensembl_ids <- paste(mitogenes$ensembl_gene_id, collapse = "|")

# Check for mitochondrial genes
# For the first SingleCellExperiment object
is_mito_raw_one <- grepl(mitogenes_ensembl_ids, rownames(sce_raw_one))
sce_raw_one$qc_stats <- perCellQCMetrics(sce_raw_one, subsets = list(Mt = is_mito_raw_one))

# For the second SingleCellExperiment object
is_mito_raw_two <- grepl(mitogenes_ensembl_ids, rownames(sce_raw_two))
sce_raw_two$qc_stats <- perCellQCMetrics(sce_raw_two, subsets = list(Mt = is_mito_raw_two))

# Eliminate samples that do not meet QC criteria
# For the first SingleCellExperiment object
libsize_drop_raw_one <- isOutlier(sce_raw_one$qc_stats$sum, nmads = 3, type = "lower", log = TRUE)
feature_drop_raw_one <- isOutlier(sce_raw_one$qc_stats$detected, nmads = 3, type = "lower", log = TRUE)
mito_drop_raw_one <- isOutlier(sce_raw_one$qc_stats$subsets_Mt_percent, nmads = 3, type = "higher")

sce_fil_raw_one <- sce_raw_one[, !(libsize_drop_raw_one | feature_drop_raw_one | mito_drop_raw_one)]

filtered_cells_raw_one <- colnames(sce_fil_raw_one)

# For the second SingleCellExperiment object
libsize_drop_raw_two <- isOutlier(sce_raw_two$qc_stats$sum, nmads = 3, type = "lower", log = TRUE)
feature_drop_raw_two <- isOutlier(sce_raw_two$qc_stats$detected, nmads = 3, type = "lower", log = TRUE)
mito_drop_raw_two <- isOutlier(sce_raw_two$qc_stats$subsets_Mt_percent, nmads = 3, type = "higher")

sce_fil_raw_two <- sce_raw_two[, !(libsize_drop_raw_two | feature_drop_raw_two | mito_drop_raw_two)]

filtered_cells_raw_two <- colnames(sce_fil_raw_two)

# Visualize different cell types in the datasets
print(unique(colData(sce_raw_one)$development_stage))
print(unique(colData(sce_raw_two)$cell_type))

# Create the SCE objects to be compared
# For sce_raw_one, filtering for "embryonic day 7"
sce_one <- sce_raw_one[, colData(sce_raw_one)$development_stage == "embryonic day 7"]

# For sce_raw_two, filtering for "undifferentiated human embryonic stem cells"
sce_two <- sce_raw_two[, colData(sce_raw_two)$cell_type == "undifferentiated human embryonic stem cells"]

# Normalize using deconvolution normalization
# For first data set
clusters_one <- quickCluster(sce_one)
sce_one <- computeSumFactors(sce_one, clusters = clusters_one)
sce_one <- logNormCounts(sce_one)

# For second data set
clusters_two <- quickCluster(sce_two)
sce_two <- computeSumFactors(sce_two, clusters = clusters_two)
sce_two <- logNormCounts(sce_two)

# Map the ensembl IDs to external gene symbols, entrez ID and gene biotype
# Use the Ensembl mart for human
ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes_human <- listAttributes(ensembl_human)

# Map the ensembl IDs of the first dataset
ensembl_id_one <- rownames(sce_one)

# Get the mapping information for the first dataset
mapping_info_one <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "gene_biotype",
    "description"
  ),
  filters = "ensembl_gene_id",
  values = ensembl_id_one,
  mart = ensembl_human
)

map_df_one <- as.data.frame(mapping_info_one)

# Plot the distribution of gene biotypes
one_gb <- ggplot(map_df_one, aes(x = gene_biotype)) +
  geom_bar(fill = "blue") +
  theme_classic() +
  labs(
    title = "Count of Genes by Gene Biotype : Dataset One",
    x = "Gene Biotype",
    y = "Count"
  ) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 10, face = "bold", hjust = 0.5, vjust = -0.5),
    axis.title.y = element_text(size = 10, face = "bold", hjust = 0.5, vjust = 1.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10)
  )

print(one_gb)

# Second
# Map the ensembl IDs of the second dataset
ensembl_id_two <- rownames(sce_two)

# Get the mapping information for the second dataset
mapping_info_two <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "gene_biotype",
    "description"
  ),
  filters = "ensembl_gene_id",
  values = ensembl_id_two,
  mart = ensembl_human
)

map_df_two <- as.data.frame(mapping_info_two)

# Plot the distribution of gene biotypes
two_gb <- ggplot(map_df_two, aes(x = gene_biotype)) +
  geom_bar(fill = "blue") +
  theme_classic() +
  labs(
    title = "Count of Genes by Gene Biotype : Dataset Two",
    x = "Gene Biotype",
    y = "Count"
  ) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 10, face = "bold", hjust = 0.5, vjust = -0.5),
    axis.title.y = element_text(size = 10, face = "bold", hjust = 0.5, vjust = 1.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10)
  )

print(two_gb)

# Identify stably expressed genes in the SCE objects
# Order the row names
sce_one <- sce_one[order(rownames(sce_one), decreasing = FALSE), ]
sce_two <- sce_two[order(rownames(sce_two), decreasing = FALSE), ]

# Identify the stably expressed genes
sce_one_seg <- scSEGIndex(exprs(sce_one), return_all = TRUE)
sce_two_seg <- scSEGIndex(exprs(sce_two), return_all = TRUE)

# Find the top 1000 stable genes in the two datasets
sce_one_seg_1000 <- sce_one_seg %>%
  arrange(desc(segIdx)) %>%
  head(1000)

sce_two_seg_1000 <- sce_two_seg %>%
  arrange(desc(segIdx)) %>%
  head(1000)

# Extract the most stable gene names from the two data sets
first_1000_genes <- rownames(as.data.frame(sce_one_seg_1000))
second_1000_genes <- rownames(as.data.frame(sce_two_seg_1000))

# Find the common genes between the two data sets
final_stable_gene_list <- union(first_1000_genes, second_1000_genes)
unique_genes <- Filter(function(x) length(grep(x, final_stable_gene_list)) == 1, final_stable_gene_list)

# Extract the stability index scores for those cumulative top stable genes
# Create data frame from the results of scSEGIndex() 
sce_one_seg_df <- as.data.frame(sce_one_seg)
sce_two_seg_df <- as.data.frame(sce_two_seg)

# Arrange the SEG tables in descending order of segIdx values
sce_one_seg_df <- sce_one_seg_df %>%
  arrange(desc(segIdx))

sce_two_seg_df <- sce_two_seg_df %>%
  arrange(desc(segIdx))

# Save the original SEG tables
write.csv(sce_one_seg_df, file = "./output_files/original_EMTAB3929_SEG_table.csv", row.names = FALSE)
write.csv(sce_two_seg_df, file = "./output_files/original_GSE75748_SEG_table.csv", row.names = FALSE)

# Only retain the rows corresponding to the most stable genes
sce_one_seg_df <- sce_one_seg_df[unique_genes,]
sce_two_seg_df <- sce_two_seg_df[unique_genes,]

# Rename the segIdx columns
names(sce_one_seg_df)[names(sce_one_seg_df) == "segIdx"] <- "segIdx_one"
names(sce_two_seg_df)[names(sce_two_seg_df) == "segIdx"] <- "segIdx_two"

# Arrange the data frame in ascending order of segIdx values
sce_one_seg_df <- sce_one_seg_df %>%
  arrange(desc(segIdx_one))

sce_two_seg_df <- sce_two_seg_df %>%
  arrange(desc(segIdx_two))

# Cross add the segIdx columns to the data frames
sce_one_seg_df$segIdx_two <- sce_two_seg_df[rownames(sce_one_seg_df), "segIdx_two"]
sce_two_seg_df$segIdx_one <- sce_one_seg_df[rownames(sce_two_seg_df), "segIdx_one"]

# Compute the differences between the segIdx columns
sce_one_seg_df$segIdx_one_vs_two <- sce_one_seg_df$segIdx_one - sce_one_seg_df$segIdx_two
sce_two_seg_df$segIdx_two_vs_one <- sce_two_seg_df$segIdx_two - sce_two_seg_df$segIdx_one

# Store the columns with NA values (NA values arose due to high percentage of zero number)
na_rows_one <- sce_one_seg_df %>%
  filter(is.na(segIdx_one_vs_two)) %>%
  rownames()

na_rows_two <- sce_two_seg_df %>%
  filter(is.na(segIdx_two_vs_one)) %>%
  rownames()

# Remove NA values from the segIdx difference columns
# Remove rows with NA values in the segIdx column
sce_one_seg_df <- sce_one_seg_df %>%
  filter(!is.na(segIdx_one_vs_two))

sce_two_seg_df <- sce_two_seg_df %>%
  filter(!is.na(segIdx_two_vs_one))

# Rename the data frames
sce_one_df <- sce_one_seg_df
sce_two_df <- sce_two_seg_df

# Add the external gene symbols to the seg data frames
# Create new data frame containing the mapping of ensembl IDs to external gene name for the first dataset
ensb_to_name_one <- data.frame(
  gene = map_df_one$ensembl_gene_id,
  external_gene_name = map_df_one$external_gene_name
)

# Create new data frame containing the mapping of ensembl IDs to external gene name for the second dataset
ensb_to_name_two <- data.frame(
  gene = map_df_two$ensembl_gene_id,
  external_gene_name = map_df_two$external_gene_name
)

# Join the data frames using an inner join
sce_one_df <- inner_join(sce_one_df, ensb_to_name_one, by = c("gene" = "gene"))
sce_two_df <- inner_join(sce_two_df, ensb_to_name_two, by = c("gene" = "gene"))

# Remove rows that don't have a valid gene name (most likely long non coding RNA genes and pseudo genes)
sce_one_df_clean <- sce_one_df %>%
  filter(!is.na(external_gene_name) & external_gene_name != "")

sce_two_df_clean <- sce_two_df %>%
  filter(!is.na(external_gene_name) & external_gene_name != "")

# Prepare the input list for GSEA
# Arrange the data frames in ascending order
sce_one_df_clean <- sce_one_df_clean %>% arrange(desc(sce_one_df_clean))
sce_two_df_clean <- sce_two_df_clean %>% arrange(desc(sce_two_df_clean))

# Save the seg tables as CSV files
write.csv(sce_one_df_clean, file = "./output_files/unionlist_EMTAB3929_SEG_table.csv", row.names = FALSE)
write.csv(sce_two_df_clean, file = "./output_files/unionlist_GSE75748_SEG_table.csv", row.names = FALSE)

# Create named vector list
sce_one_lst <- setNames(sce_one_df_clean$segIdx_one_vs_two, sce_one_df_clean$external_gene_name)
sce_two_lst <- setNames(sce_two_df_clean$segIdx_two_vs_one, sce_two_df_clean$external_gene_name)

# Retrieve the pathway lists from molecular signatures database
# Retrieve the MSigDB (Molecular Signatures Database) for Homo sapiens (human) using gene symbols (version 7.5.1)
# SYM is for gene symbols other possible values: 'ENTREZ', 'ENSEMBL', 'REFSEQ', 'UNIPROT', 'HGNC'
msigdb.hs = getMsigdb(org = 'hs', id = 'SYM', version = '7.5.1')
msigdb.hs = appendKEGG(msigdb.hs)
msigdb.hs.list <- geneIds(msigdb.hs)

save(msigdb.hs.list, file = "./final_msigdb_hs_list.RData")

# Perform GSEA using fgsea for all pathways
# for the first subset
fgsea_one <- fgsea(pathways = msigdb.hs.list, stats = sce_one_lst, minSize = 5, maxSize = 1000, nPermSimple = 100000)
fgsea_one <- fgsea_one %>%
  arrange(padj, desc(abs(NES)))

# for the second subset
fgsea_two <- fgsea(pathways = msigdb.hs.list, stats = sce_two_lst, minSize = 5, maxSize = 1000, nPermSimple = 100000)
fgsea_two <- fgsea_two %>%
  arrange(padj, desc(abs(NES)))

# Perform GSEA using fgsea for hallmark pathways (optional)
# Filter for Hallmark gene sets
hallmark_gene_sets_list <- msigdb.hs.list[grep("^HALLMARK_", names(msigdb.hs.list))]

# for the first subset
hallmark_fgsea_one <- fgsea(pathways = hallmark_gene_sets_list, stats = sce_one_lst, minSize = 5, maxSize = 1000, nPermSimple = 100000)
hallmark_fgsea_one <- hallmark_fgsea_one %>%
  arrange(padj, desc(abs(NES)))

# for the second subset
hallmark_fgsea_two <- fgsea(pathways = hallmark_gene_sets_list, stats = sce_two_lst, minSize = 5, maxSize = 1000, nPermSimple = 100000)
hallmark_fgsea_two <- hallmark_fgsea_two %>%
  arrange(padj, desc(abs(NES)))

# Save the GSEA results
# Function to convert list columns to JSON strings
convert_list_columns_to_json <- function(df) {
  list_columns <- sapply(df, is.list)
  df[list_columns] <- lapply(df[list_columns], function(col) {
    sapply(col, jsonlite::toJSON, auto_unbox = TRUE)
  })
  return(df)
}

# Convert the list to dataframe
fgsea_one <- convert_list_columns_to_json(as.data.frame(fgsea_one))
fgsea_two <- convert_list_columns_to_json(as.data.frame(fgsea_two))
hallmark_fgsea_one <- convert_list_columns_to_json(as.data.frame(hallmark_fgsea_one))
hallmark_fgsea_two <- convert_list_columns_to_json(as.data.frame(hallmark_fgsea_two))

# Save as csv files
write.csv(fgsea_one, file = "./output_files/EMTAB3929vsGSE75748_fgsea.csv", row.names = FALSE)
write.csv(fgsea_two, file = "./output_files/GSE75748vsEMTAB3929_fgsea.csv", row.names = FALSE)
write.csv(hallmark_fgsea_one, file = "./output_files/EMTAB3929vsGSE75748_fgsea_hallmark.csv", row.names = FALSE)
write.csv(hallmark_fgsea_two, file = "./output_files/GSE75748vsEMTAB3929_fgsea_hallmark.csv", row.names = FALSE)

# Downstream analysis

# Plot enrichment for first pathway
plotEnrichment(msigdb.hs.list[["BENPORATH_ES_1"]], sce_one_lst) +
  ggtitle(paste("Enrichment plot for BENPORATH_ES_1 in E-MTAB-3929"))

plotEnrichment(msigdb.hs.list[["BENPORATH_ES_1"]], sce_two_lst) +
  ggtitle(paste("Enrichment plot for BENPORATH_ES_1 in GSE75748"))

# Combine plots for visualization
library(cowplot)

plot1 <- plotEnrichment(msigdb.hs.list[["BENPORATH_ES_1"]], sce_one_lst) +
  ggtitle("Enrichment plot for BENPORATH_ES_1 in E-MTAB-3929")

plot2 <- plotEnrichment(msigdb.hs.list[["BENPORATH_ES_1"]], sce_two_lst) +
  ggtitle("Enrichment plot for BENPORATH_ES_1 in GSE75748")

combined_plot <- plot_grid(plot1, plot2, labels = c("A", "B"), ncol = 2)

print(combined_plot)

# For Second pathway
plot3 <- plotEnrichment(msigdb.hs.list[["BENPORATH_NANOG_TARGETS"]], sce_one_lst) +
  ggtitle("Enrichment plot for BENPORATH_NANOG_TARGETS in E-MTAB-3929")

plot4 <- plotEnrichment(msigdb.hs.list[["BENPORATH_NANOG_TARGETS"]], sce_two_lst) +
  ggtitle("Enrichment plot for BENPORATH_NANOG_TARGETS in GSE75748")

combined_plot_1 <- plot_grid(plot3, plot4, labels = c("C", "D"), ncol = 2)

print(combined_plot_1)

# For OCT4 targets
plot5 <- plotEnrichment(msigdb.hs.list[["BENPORATH_OCT4_TARGETS"]], sce_one_lst) +
  ggtitle("Enrichment plot for BENPORATH_OCT4_TARGETS in E-MTAB-3929")

plot6 <- plotEnrichment(msigdb.hs.list[["BENPORATH_OCT4_TARGETS"]], sce_two_lst) +
  ggtitle("Enrichment plot for BENPORATH_OCT4_TARGETS in GSE75748")

combined_plot_2 <- plot_grid(plot5, plot6, labels = c("A", "B"), ncol = 2)

print(combined_plot_2)

# For SOX2 targets
plot7 <- plotEnrichment(msigdb.hs.list[["BENPORATH_SOX2_TARGETS"]], sce_one_lst) +
  ggtitle("Enrichment plot for BENPORATH_SOX2_TARGETS in E-MTAB-3929")

plot8 <- plotEnrichment(msigdb.hs.list[["BENPORATH_SOX2_TARGETS"]], sce_two_lst) +
  ggtitle("Enrichment plot for BENPORATH_SOX2_TARGETS in GSE75748")

combined_plot_3 <- plot_grid(plot7, plot8, labels = c("E", "F"), ncol = 2)

print(combined_plot_3)

# Final combined plot for OCT4, NANOG and SOX2 targets
final_combined_plot <- plot_grid(combined_plot_2, combined_plot_1, combined_plot_3, nrow = 3)
ggsave("./figures/combined_enrichment_plot.png", final_combined_plot, width = 16, height = 18, dpi = 500)
print(final_combined_plot)

# Save session data and history
save.image(file = "finalized-human-ranked-SEG-enrichment.RData")
savehistory(file = "finalized-human-ranked-SEG-enrichment")
writeLines(capture.output(sessionInfo()), "finalized-human-ranked-SEG-enrichment_session_info.txt")

# Check the position of SOX2, NANOG and Prdm14 in the output dataframe of scSEGIndex()
sox2 <- "ENSG00000181449"
nanog <- "ENSG00000111704"
prdm14 <- "ENSG00000147596"

genes_to_check <- c(sox2, nanog, prdm14)

# Arrange the data frames by segIdx
seg_one_check_df <- sce_one_seg %>% arrange(segIdx)
seg_two_check_df <- sce_two_seg %>% arrange(segIdx)

# Find the row indices where the values are present in each data frame and set the names
row_indices_sce_one <- setNames(which(seg_one_check_df$gene %in% genes_to_check), seg_one_check_df$gene[which(seg_one_check_df$gene %in% genes_to_check)])
row_indices_sce_two <- setNames(which(seg_two_check_df$gene %in% genes_to_check), seg_two_check_df$gene[which(seg_two_check_df$gene %in% genes_to_check)])

print(row_indices_sce_one)
print(row_indices_sce_two)

# Visualize the expression of these genes in the datasets
# Combined heatmaps
# Define function for heat map creation
create_heatmap <- function(sce, gene, main_title = "") {
  normalized_counts <- assay(sce, "logcounts")
  
  if (!(gene %in% rownames(normalized_counts))) {
    stop(paste("Gene", gene, "not found in the dataset."))
  }
  
  gene_expression <- normalized_counts[gene, , drop = FALSE]
  gene_expression <- t(gene_expression)
  gene_expression <- scale(gene_expression, center = TRUE, scale = TRUE)
  
  pheatmap(gene_expression,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = FALSE,
           show_colnames = TRUE,
           main = main_title,
           silent = TRUE)
}

# Define genes for which the heatmap comparisons are to be made
# Define the genes
gene1 <- "ENSG00000111704" # NANOG
gene2 <- "ENSG00000147596" # PRDM14
gene3 <- "ENSG00000181449" # SOX2
gene4 <- "ENSG00000104332" # SFRP1
gene5 <- "ENSG00000152284" # TCF7L1
gene6 <- "ENSG00000204531" # POU5F1 (has got multiple ensembl IDs)

# NANOG
# Generate heat maps
heatmap1 <- create_heatmap(sce_one, gene1)
heatmap2 <- create_heatmap(sce_two, gene1)

# Combine heat maps using cow plot
nanog_plot <- plot_grid(heatmap1$gtable, heatmap2$gtable, labels = c("A", "B"))

# Display the combined plot
print(nanog_plot)

# PRDM14
# Generate heat maps
heatmap3 <- create_heatmap(sce_one, gene2)
heatmap4 <- create_heatmap(sce_two, gene2)

# Combine heat maps using cow plot
prdm14_plot <- plot_grid(heatmap3$gtable, heatmap4$gtable, labels = c("A", "B"))

print(prdm14_plot)

# SOX2
# Generate heat maps
heatmap5 <- create_heatmap(sce_one, gene3)
heatmap6 <- create_heatmap(sce_two, gene3)

# Combine heat maps using cow plot
sox2_plot <- plot_grid(heatmap5$gtable, heatmap6$gtable, labels = c("A", "B"))

print(sox2_plot)

# SFRP1, TCF7L1 and POU5F1
# Generate heat maps
heatmap7 <- create_heatmap(sce_one, gene4, main_title = "Expression of SFRP1 in E-MTAB-3929")
heatmap8 <- create_heatmap(sce_two, gene4, main_title = "Expression of SFRP1 in GSE75748")

heatmap9 <- create_heatmap(sce_one, gene5, main_title = "Expression of TCF5L1 in E-MTAB-3929")
heatmap10 <- create_heatmap(sce_two, gene5, main_title = "Expression of TCF5L1 in GSE75748")

heatmap11 <- create_heatmap(sce_one, gene6, main_title = "Expression of POU5F1 in E-MTAB-3929")
heatmap12 <- create_heatmap(sce_two, gene6, main_title = "Expression of POU5F1 in GSE75748")

# Combine heat maps using cow plot
sfrp1_plot <- plot_grid(heatmap7$gtable, heatmap8$gtable, labels = c("A", "B"))
tcf7l1_plot <- plot_grid(heatmap9$gtable, heatmap10$gtable, labels = c("A", "B"))
pou5f1_plot <- plot_grid(heatmap11$gtable, heatmap12$gtable, labels = c("A", "B"))

# Define a function to save plots
save_combined_plot <- function(plot, filename) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  ggsave(filename, plot, width = 10, height = 10, units = "in", dpi = 500)
}

save_combined_plot(sfrp1_plot, "./figures/sfrp1_heatmap.png")
save_combined_plot(tcf7l1_plot, "./figures/tcf7l1_heatmap.png")
save_combined_plot(pou5f1_plot, "./figures/pou5f1_heatmap.png")

# Heatmaps for most stable and least stable genes
# Identify the top 20 most stable and least stable protein coding genes
# Arrange the original SEG table in descending order of segIdx column
seg_table_one <- as.data.frame(sce_one_seg) %>% arrange(desc(segIdx))
seg_table_two <- as.data.frame(sce_two_seg) %>% arrange(desc(segIdx))

# Map the row names of the seg table to external gene names and gene bio types
seg_table_one <- inner_join(seg_table_one, map_df_one, by = c("gene" = "ensembl_gene_id"))
seg_table_two <- inner_join(seg_table_two, map_df_two, by = c("gene" = "ensembl_gene_id"))

# Subset the table to only contain protein coding genes
seg_table_one_pc <- seg_table_one %>% filter(gene_biotype == "protein_coding")
seg_table_two_pc <- seg_table_two %>% filter(gene_biotype == "protein_coding")

# Arrange the new data frame in descending order of the segIdx column
seg_table_one_pc <- seg_table_one_pc %>% arrange(desc(segIdx))
seg_table_two_pc <- seg_table_two_pc %>% arrange(desc(segIdx))

# Convert the objects to data frames and drop rows with NA in segIdx column
seg_table_one_pc_df <- as.data.frame(seg_table_one_pc) %>% drop_na(segIdx)
seg_table_two_pc_df <- as.data.frame(seg_table_two_pc) %>% drop_na(segIdx)

# Extract the top 20 most stable and least stable rows based on the gene column
top_20_most_stable_one <- head(seg_table_one_pc_df, 20)$gene
top_20_least_stable_one <- tail(seg_table_one_pc_df, 20)$gene
top_20_most_stable_two <- head(seg_table_two_pc_df, 20)$gene
top_20_least_stable_two <- tail(seg_table_two_pc_df, 20)$gene

# Define the heat map function
create_heatmap <- function(expression_matrix, genes, map_df_one, sce_title = NULL) {
  # Check if all genes are in the dataset
  missing_genes <- setdiff(genes, rownames(expression_matrix))
  if (length(missing_genes) > 0) {
    stop(paste("Genes not found in the dataset:", paste(missing_genes, collapse = ", ")))
  }
  
  # Extract gene expression for selected genes
  gene_expression <- expression_matrix[genes, , drop = FALSE]
  
  # Column normalization: center and scale each column to have mean = 0 and sd = 1
  gene_expression <- t(scale(t(gene_expression), center = TRUE, scale = TRUE))
  
  # Map Ensembl IDs to external gene names using map_df_one
  external_gene_names <- map_df_one$external_gene_name[match(genes, map_df_one$ensembl_gene_id)]
  
  # Identify and report any Ensembl IDs that could not be mapped
  if (any(is.na(external_gene_names))) {
    unmapped_genes <- genes[is.na(external_gene_names)]
    stop(paste("Some Ensembl IDs could not be mapped to external gene names:", paste(unmapped_genes, collapse = ", ")))
  }
  
  # Transpose for heatmap visualization
  gene_expression <- t(gene_expression)
  
  # Set the column names to external gene names
  colnames(gene_expression) <- external_gene_names
  
  # Set the title
  title <- if (!is.null(sce_title)) {
    paste("Expression of selected genes in", sce_title)
  } else {
    "Expression of selected genes"
  }
  
  # Generate the heatmap
  heatmap <- pheatmap(gene_expression,
                      cluster_rows = TRUE,
                      cluster_cols = FALSE,
                      show_rownames = FALSE,
                      show_colnames = TRUE,
                      main = title)
  
  return(heatmap)
}

# Extract the normalized counts for the data sets
expression_matrix_one <- assay(sce_one, "logcounts")
expression_matrix_two <- assay(sce_two, "logcounts")

# Plot the heatmaps for visualizing the expression patterns of top 20 most and least stable genes in the data sets
# Top 20 most stable genes in E-MTAB-3929
create_heatmap(expression_matrix_one, top_20_most_stable_one, sce_title = "E-MTAB-3929 - Most Stable")

# Top 20 least stable genes in E-MTAB-3929
create_heatmap(expression_matrix_one, top_20_least_stable_one, sce_title = "E-MTAB-3929 - Least Stable")

# Top 20 most stable genes in GSE75748
create_heatmap(expression_matrix_two, top_20_most_stable_two, sce_title = "GSE75748 - Most Stable")

# Top 20 least stable genes in GSE75748
create_heatmap(expression_matrix_two, top_20_least_stable_two, sce_title = "GSE75748 - Least Stable")

# Save all the heatmaps as high quality PNG file
# Define function for saving the heat maps
save_pheatmap_png <- function(heatmap, filename) {
  png(filename, width = 8, height = 8, units = "in", res = 500)
  grid::grid.draw(heatmap$gtable)
  dev.off()
}

# Create heat maps and store them in variables
heatmap_most_stable_one <- create_heatmap(expression_matrix_one, top_20_most_stable_one, map_df_one, sce_title = "E-MTAB-3929 - Most Stable")
heatmap_least_stable_one <- create_heatmap(expression_matrix_one, top_20_least_stable_one, map_df_one, sce_title = "E-MTAB-3929 - Least Stable")
heatmap_most_stable_two <- create_heatmap(expression_matrix_two, top_20_most_stable_two, map_df_one, sce_title = "GSE75748 - Most Stable")
heatmap_least_stable_two <- create_heatmap(expression_matrix_two, top_20_least_stable_two, map_df_one, sce_title = "GSE75748 - Least Stable")

save_pheatmap_png(heatmap_most_stable_one, "./Heatmaps/E-MTAB-3929_Most_Stable.png")
save_pheatmap_png(heatmap_least_stable_one, "./Heatmaps/E-MTAB-3929_Least_Stable.png")
save_pheatmap_png(heatmap_most_stable_two, "./Heatmaps/GSE75748_Most_Stable.png")
save_pheatmap_png(heatmap_least_stable_two, "./Heatmaps/GSE75748_Least_Stable.png")
