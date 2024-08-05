# Load the human stably expressed gene from scMerge package
data("segList_ensemblGeneID", package = "scMerge")
hSEG_original <- segList_ensemblGeneID$human$human_scSEG

# Create a data frame containing the rank
hSEG_original_df <- tibble(gene = hSEG_original) %>%
  mutate(rank = row_number())

# Use biomaRt to map Ensembl IDs to external gene symbols and gene bio types
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve the external gene symbols and gene bio types
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = hSEG_original_df$gene,
  mart = mart
)

# Merge the retrieved information with the original data frame
hSEG_original_df <- hSEG_original_df %>%
  left_join(gene_info, by = c("gene" = "ensembl_gene_id"))

## Now find the position of these ensembl Ids in the seg tables computed by scSEGIndex() in both the data sets

# Arrange the seg table in descending order of segIdx value
first_seg_table <- sce_one_seg %>% arrange(desc(segIdx))
second_seg_table <- sce_two_seg %>% arrange(desc(segIdx))


# Find the positions of the Ensembl IDs in the first_seg_table
positions_one <- match(hSEG_original_df$gene, first_seg_table$gene)
result_df_one <- tibble(
  gene = hSEG_original_df$gene,
  position = positions_one
)

# Find the positions of the Ensembl IDs in the second_seg_table
positions_two <- match(hSEG_original_df$gene, second_seg_table$gene)
result_df_two <- tibble(
  gene = hSEG_original_df$gene,
  position = positions_two
)

# In the results data frame map the ensembl IDs to gene names
result_df_one_ordered <- result_df_one %>%
  left_join(gene_info, by = c("gene" = "ensembl_gene_id")) %>%
  arrange(position)

result_df_two_ordered <- result_df_two %>%
  left_join(gene_info, by = c("gene" = "ensembl_gene_id")) %>%
  arrange(position)

# Save the results as csv files
write.csv(result_df_one_ordered, file = "./hSEG_position_EMTAB3929.csv", row.names = FALSE)
write.csv(result_df_two_ordered, file = "./hSEG_position_GSE75748.csv", row.names = FALSE)