# Order the rownames of sce_hESC
sce_hESC <- sce_hESC[order(rownames(sce_hESC), decreasing = FALSE), ]

# Order the rownames of sce_non_hESC
sce_non_hESC <- sce_non_hESC[order(rownames(sce_non_hESC), decreasing = FALSE), ]

# Calculate Stability Indices
res_hESC <- scSEGIndex(exprs(sce_hESC), return_all = TRUE)
res_non_hESC <- scSEGIndex(exprs(sce_non_hESC), return_all = TRUE)

res_hESC_ordered <- res_hESC[order(res_hESC$segIdx, decreasing = FALSE), ]
res_non_hESC_ordered <- res_non_hESC[order(res_non_hESC$segIdx, decreasing = FALSE), ]

res_hESC_ordered$segIdx_hESC <- rank(-1 * res_hESC_ordered$segIdx, ties.method = "average")
res_non_hESC_ordered$segIdx_non_hESC <- rank(-1 * res_non_hESC_ordered$segIdx, ties.method = "average")

res_hESC_ordered$segIdx_non_hESC <- res_non_hESC_ordered[rownames(res_hESC_ordered), ]$segIdx_non_hESC
res_non_hESC_ordered$segIdx_hESC <- res_hESC_ordered[rownames(res_non_hESC_ordered), ]$segIdx_hESC

res_hESC_ordered$segIdx_hESCvnon_hESC <- res_hESC_ordered$segIdx_hESC - res_hESC_ordered$segIdx_non_hESC
res_non_hESC_ordered$segIdx_non_hESCvhESC <- res_non_hESC_ordered$segIdx_non_hESC - res_hESC_ordered$segIdx_hESC

res_hESC_ordered <- res_hESC_ordered[-12]

res_hESC_ordered <- res_hESC_ordered[order(res_hESC_ordered$segIdx, decreasing = TRUE), ]
res_non_hESC_ordered <- res_non_hESC_ordered[order(res_non_hESC_ordered$segIdx, decreasing = TRUE), ]

res_hESC_ordered1000 <- res_hESC_ordered[1:1000, ]
res_non_hESC_ordered1000 <- res_non_hESC_ordered[1:1000, ]

res_hESC_ordered1000 <- res_hESC_ordered1000[order(res_hESC_ordered1000$segIdx_hESCvnon_hESC, decreasing = TRUE), ]
res_non_hESC_ordered1000 <- res_non_hESC_ordered1000[order(res_non_hESC_ordered1000$segIdx_non_hESCvhESC, decreasing = TRUE), ]

res_hESC_ordered1000_lst <- res_hESC_ordered1000$segIdx_hESCvnon_hESC
names(res_hESC_ordered1000_lst) <- rownames(res_hESC_ordered1000)

res_non_hESC_ordered1000_lst <- res_non_hESC_ordered1000$segIdx_non_hESCvhESC
names(res_non_hESC_ordered1000_lst) <- rownames(res_non_hESC_ordered1000)

# Convert Ensembl IDs to gene symbols
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_annotations <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                          filters = 'ensembl_gene_id',
                          values = rownames(res_hESC_ordered1000),
                          mart = ensembl)

res_hESC_ordered1000$external_gene_name <- gene_annotations$external_gene_name[match(rownames(res_hESC_ordered1000), gene_annotations$ensembl_gene_id)]
res_non_hESC_ordered1000$external_gene_name <- gene_annotations$external_gene_name[match(rownames(res_non_hESC_ordered1000), gene_annotations$ensembl_gene_id)]

# Perform functional enrichment analysis
eh <- ExperimentHub()
query(eh, 'msigdb')

msigdb.hs <- getMsigdb(org = 'hs', id = 'SYM', version = '7.5.1')
msigdb.hs <- appendKEGG(msigdb.hs)
msigdb.hs.list <- geneIds(msigdb.hs)

# Ensure that the external gene names are correctly mapped
res_hESC_ordered1000$external_gene_name <- gene_annotations$external_gene_name[match(rownames(res_hESC_ordered1000), gene_annotations$ensembl_gene_id)]
res_non_hESC_ordered1000$external_gene_name <- gene_annotations$external_gene_name[match(rownames(res_non_hESC_ordered1000), gene_annotations$ensembl_gene_id)]

# Ensure no NA values remain
res_hESC_ordered1000 <- res_hESC_ordered1000[!is.na(res_hESC_ordered1000$external_gene_name), ]
res_non_hESC_ordered1000 <- res_non_hESC_ordered1000[!is.na(res_non_hESC_ordered1000$external_gene_name), ]

# Prepare the ranked lists
res_hESC_ordered1000_lst <- setNames(res_hESC_ordered1000$segIdx_hESCvnon_hESC, res_hESC_ordered1000$external_gene_name)
res_non_hESC_ordered1000_lst <- setNames(res_non_hESC_ordered1000$segIdx_non_hESCvhESC, res_non_hESC_ordered1000$external_gene_name)

# Filter Pathways using gene symbols based on size
min_size <- 5
max_size <- 500
filtered_pathways <- Filter(function(p) length(p) >= min_size & length(p) <= max_size, msigdb.hs.list)

# Verify gene overlap with pathways using gene symbols
gene_set_hESC <- names(res_hESC_ordered1000_lst)
gene_set_non_hESC <- names(res_non_hESC_ordered1000_lst)

overlap_hESC <- sapply(filtered_pathways, function(p) length(intersect(p, gene_set_hESC)))
overlap_non_hESC <- sapply(filtered_pathways, function(p) length(intersect(p, gene_set_non_hESC)))

print(table(overlap_hESC))
print(table(overlap_non_hESC))

# Perform fgsea with adjusted parameters
fgsea_res_hESC <- fgsea(filtered_pathways, res_hESC_ordered1000_lst, nPermSimple = 200000, minSize = min_size, maxSize = max_size)
fgsea_res_non_hESC <- fgsea(filtered_pathways, res_non_hESC_ordered1000_lst, nPermSimple = 200000, minSize = min_size, maxSize = max_size)

# Check results
print(fgsea_res_hESC)
print(fgsea_res_non_hESC)
fgsea_res_non_hESC[1:20,]$leadingEdge
fgsea_res_hESC[1:20,]$leadingEdge

fgsea_hESC_ordered1000_lstr <- fgsea_res_hESC[order(padj, -abs(NES)), ]
fgsea_hESC_ordered1000_lstr[1:20,]
fgsea_hESC_ordered1000_lstr[1:20,]$leadingEdge