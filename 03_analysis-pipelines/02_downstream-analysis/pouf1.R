library(biomaRt)
library(tidyverse)


# Function to retrieve Ensembl IDs for a given gene symbol
get_ensembl_ids <- function(gene_symbol) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  results <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name","gene_biotype"),
    filters = "external_gene_name",
    values = gene_symbol,
    mart = ensembl
  )
  
  return(results)
}


gene_symbol <- "POU5F1"
ensembl_ids <- get_ensembl_ids(gene_symbol)
print(ensembl_ids)

gene_symbol_df <- as.data.frame(ensembl_ids)
pou5f1_list <- gene_symbol_df$ensembl_gene_id

# Check which ensembl Ids of pou5f1 is in rownames of counts matrix and in hSEG list
present_in_counts <- pou5f1_list %in% rownames(raw_one)


data("segList_ensemblGeneID", package = "scMerge")
hSEG <- segList_ensemblGeneID$human$human_scSEG


present_in_hSEG <- pou5f1_list %in% hSEG

# No variant of pou5f1 is present in the hSEG list

# Only one ensembl ID for POU5F1 present in the counts matrix and it is as follows
print(pou5f1_list[present_in_counts])

