# Intitial analysis
library(Seurat)
library(scMerge)
library(SingleCellExperiment)
library(ggplot2)
library(tidyverse)
library(biomaRt)
library(corrplot)
library(scater)
library(scran)
library(biomaRt)
# GSEA and Functional Enrichment
library(msigdb)
library(fgsea)
library(ExperimentHub)
library(GSEABase)


# Load in expression data and annotation data
raw <- readRDS("../../study-datasets/mouse-dataset/E-MTAB-2600/single_es_tximport.rds")
annot <- read.table("../../study-datasets/mouse-dataset/E-MTAB-2600/E-MTAB-2600.targets.txt", sep = "\t", header = TRUE, row.names = "ERR")
annot <- annot[match(colnames(raw$counts), rownames(annot)), ]


# Load in SEG
data("segList_ensemblGeneID", package = "scMerge")
mseg <- segList_ensemblGeneID$mouse$mouse_scSEG

# SCE object
sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(raw$counts)),
  colData = annot
)

# Sample Quality Control
is_spike_one <- grepl("^ERCC", rownames(sce))
mitogenes <- "ENSMUSG00000064336|ENSMUSG00000064337|ENSMUSG00000064338|ENSMUSG00000064339|ENSMUSG00000064340|ENSMUSG00000064341|ENSMUSG00000064342|ENSMUSG00000064343|ENSMUSG00000064344|ENSMUSG00000064345|ENSMUSG00000064346|ENSMUSG00000064347|ENSMUSG00000064348|ENSMUSG00000064349|ENSMUSG00000064350|ENSMUSG00000064351|ENSMUSG00000064352|ENSMUSG00000064353|ENSMUSG00000064354|ENSMUSG00000064355|ENSMUSG00000064356|ENSMUSG00000064357|ENSMUSG00000064358|ENSMUSG00000064359|ENSMUSG00000064360|ENSMUSG00000064361|ENSMUSG00000065947|ENSMUSG00000064363|ENSMUSG00000064364|ENSMUSG00000064365|ENSMUSG00000064366|ENSMUSG00000064367|ENSMUSG00000064368|ENSMUSG00000064369|ENSMUSG00000064370|ENSMUSG00000064371|ENSMUSG00000064372"
is_mito_one <- grepl(mitogenes, rownames(sce))
sce$qc_stats <- perCellQCMetrics(sce, subsets=list(ERCC=is_spike_one, Mt=is_mito_one))

# Calculate the QC metrics based on which cells will be retained or eliminated
libsize_drop_one <- isOutlier(sce$qc_stats$sum, nmads = 3, type = "lower", log = TRUE)
feature_drop_one <- isOutlier(sce$qc_stats$detected, nmads = 3, type = "lower", log = TRUE)
mito_drop_one <- isOutlier(sce$qc_stats$subsets_Mt_percent, nmads = 3, type = "higher")
spike_drop_one <- isOutlier(sce$qc_stats$subsets_ERCC_percent, nmads = 3, type = "higher")

# Save a copy of the unfiltered SCE object
sce_copy <- sce

# Drop the samples from the SCE object that do not meet the QC standards
sce_fil <- sce[, !(libsize_drop_one | feature_drop_one | mito_drop_one | spike_drop_one)]
filtered_cells_one <- colnames(sce_fil)

# Subset the SCE object for 2i and serum cell types
sce_2i <- sce_fil[, colData(sce_fil)$Type == "2i"]
sce_serum <- sce_fil[, colData(sce_fil)$Type == "serum"]

# Normalize both the SCE objects using deconvolution normalization
# For 2i subset
clusters_2i <- quickCluster(sce_2i)
sce_2i <- computeSumFactors(sce_2i, clusters = clusters_2i)
sce_2i <- logNormCounts(sce_2i)

# For serum subset
clusters_serum <- quickCluster(sce_serum)
sce_serum <- computeSumFactors(sce_serum, clusters = clusters_serum)
sce_serum <- logNormCounts(sce_serum)

# Order the rownames of sce_serum
sce_serum <- sce_serum[order(rownames(sce_serum), decreasing = FALSE), ]
# Order the rownames of sce_2i
sce_2i <- sce_2i[order(rownames(sce_2i), decreasing = FALSE), ]

# Rename the variables
sce_escounts_2i <- sce_2i
sce_escounts_serum <- sce_serum




# Calculate Stability Indices
res_serum <- scSEGIndex(exprs(sce_escounts_serum), return_all = TRUE)
res_2i <- scSEGIndex(exprs(sce_escounts_2i), return_all = TRUE)

res_serum_ordered <- res_serum[order(res_serum$segIdx, decreasing = FALSE), ]
res_2i_ordered <- res_2i[order(res_2i$segIdx, decreasing = FALSE), ]

res_serum_ordered$segIdx_serum <- rank(-1 * res_serum_ordered$segIdx, ties.method = "average")
res_2i_ordered$segIdx_2i <- rank(-1 * res_2i_ordered$segIdx, ties.method = "average")

res_serum_ordered$segIdx_2i <- res_2i_ordered[rownames(res_serum_ordered), ]$segIdx_2i
res_2i_ordered$segIdx_serum <- res_serum_ordered[rownames(res_2i_ordered), ]$segIdx_serum

res_serum_ordered$segIdx_serumv2i <- res_serum_ordered$segIdx_serum - res_serum_ordered$segIdx_2i
res_2i_ordered$segIdx_2ivserum <- res_2i_ordered$segIdx_2i - res_serum_ordered$segIdx_serum

res_serum_ordered <- res_serum_ordered[-12]

res_2i_ordered <- res_2i_ordered[order(res_2i_ordered$segIdx, decreasing = TRUE), ]
res_serum_ordered <- res_serum_ordered[order(res_serum_ordered$segIdx, decreasing = TRUE), ]

res_2i_ordered1000 <- res_2i_ordered[1:1000, ]
res_serum_ordered1000 <- res_serum_ordered[1:1000, ]

res_2i_ordered1000 <- res_2i_ordered1000[order(res_2i_ordered1000$segIdx_2ivserum, decreasing = TRUE), ]
res_serum_ordered1000 <- res_serum_ordered1000[order(res_serum_ordered1000$segIdx_serumv2i, decreasing = TRUE), ]

res_2i_ordered1000_lst <- res_2i_ordered1000$segIdx_2ivserum
names(res_2i_ordered1000_lst) <- rownames(res_2i_ordered1000)

res_serum_ordered1000_lst <- res_serum_ordered1000$segIdx_serumv2i
names(res_serum_ordered1000_lst) <- rownames(res_serum_ordered1000)

# Convert Ensembl IDs to gene symbols
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
gene_annotations <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                          filters = 'ensembl_gene_id',
                          values = rownames(res_2i_ordered1000),
                          mart = ensembl)

res_2i_ordered1000$external_gene_name <- gene_annotations$external_gene_name[match(rownames(res_2i_ordered1000), gene_annotations$ensembl_gene_id)]
res_serum_ordered1000$external_gene_name <- gene_annotations$external_gene_name[match(rownames(res_serum_ordered1000), gene_annotations$ensembl_gene_id)]

#Perform functional enrichment analysis
eh <- ExperimentHub()
query(eh, 'msigdb')

msigdb.mm <- getMsigdb(org = 'mm', id = 'SYM', version = '7.5.1')
msigdb.mm <- appendKEGG(msigdb.mm)
msigdb.mm.list <- geneIds(msigdb.mm)

# Ensure that the external gene names are correctly mapped
res_2i_ordered1000$external_gene_name <- gene_annotations$external_gene_name[match(rownames(res_2i_ordered1000), gene_annotations$ensembl_gene_id)]
res_serum_ordered1000$external_gene_name <- gene_annotations$external_gene_name[match(rownames(res_serum_ordered1000), gene_annotations$ensembl_gene_id)]

#Ensure no NA values remain
res_2i_ordered1000 <- res_2i_ordered1000[!is.na(res_2i_ordered1000$external_gene_name), ]
res_serum_ordered1000 <- res_serum_ordered1000[!is.na(res_serum_ordered1000$external_gene_name), ]

#Prepare the ranked lists
res_2i_ordered1000_lst <- setNames(res_2i_ordered1000$segIdx_2ivserum, res_2i_ordered1000$external_gene_name)
res_serum_ordered1000_lst <- setNames(res_serum_ordered1000$segIdx_serumv2i, res_serum_ordered1000$external_gene_name)

#Filter Pathways using gene symbols based on size
min_size <- 15
max_size <- 500
filtered_pathways <- Filter(function(p) length(p) >= min_size & length(p) <= max_size, msigdb.mm.list)

#Verify gene overlap with pathways using gene symbols
gene_set_2i <- names(res_2i_ordered1000_lst)
gene_set_serum <- names(res_serum_ordered1000_lst)

overlap_2i <- sapply(filtered_pathways, function(p) length(intersect(p, gene_set_2i)))
overlap_serum <- sapply(filtered_pathways, function(p) length(intersect(p, gene_set_serum)))

print(table(overlap_2i))
print(table(overlap_serum))

#Perform fgsea with adjusted parameters
fgsea_res_2i <- fgsea(filtered_pathways, res_2i_ordered1000_lst, nPermSimple = 200000, minSize = min_size, maxSize = max_size)
fgsea_res_serum <- fgsea(filtered_pathways, res_serum_ordered1000_lst, nPermSimple = 200000, minSize = min_size, maxSize = max_size)

#Check results
print(fgsea_res_2i)
print(fgsea_res_serum)
fgsea_res_serum[1:20,]$leadingEdge
fgsea_res_2i[1:20,]$leadingEdge

fgsea_serum_ordered1000_lstr <- fgsea_res_serum[order(padj, -abs(NES)), ]
fgsea_serum_ordered1000_lstr[1:20,]
fgsea_serum_ordered1000_lstr[1:20,]$leadingEdge

#Save to CSV
fgsea_res_2i_dt <- as.data.table(fgsea_res_2i)
fgsea_res_2i_dt[, leadingEdge := sapply(leadingEdge, toString)]
fgsea_res_serum_dt <- as.data.table(fgsea_res_serum)
fgsea_res_serum_dt[, leadingEdge := sapply(leadingEdge, toString)]

write.csv(fgsea_res_2i_dt, "FE_results_2i.csv", row.names = FALSE)
write.csv(fgsea_res_serum_dt, "FE_results_serum.csv", row.names = FALSE)
