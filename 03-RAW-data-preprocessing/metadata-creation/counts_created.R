library(dplyr)

# Define a function to read a single file and extract relevant data
read_expression_file <- function(filepath) {
  # Read the file without specifying column names
  data <- read.table(filepath, header = FALSE, sep = "\t", fill = TRUE)
  # Select the first two columns and rename them
  gene_counts <- data.frame(gene_name = data[, 1], counts = data[, 2])
  return(gene_counts)
}

# List of all file paths
file_paths <- c(
  "GSM896803_Sample_M21.expression.txt", "GSM922195_Sample_G3.expression.txt",
  "GSM896804_Sample_M22.expression.txt", "GSM922196_Sample_G4.expression.txt",
  "GSM896805_Sample_M23.expression.txt", "GSM922197_Sample_G5.expression.txt",
  "GSM896806_Sample_C1.expression.txt", "GSM922198_Sample_G6.expression.txt",
  "GSM896807_Sample_C2.expression.txt", "GSM922199_Sample_G7.expression.txt",
  "GSM896808_Sample_C3.expression.txt", "GSM922200_Sample_G8.expression.txt",
  "GSM896809_Sample_B1.expression.txt", "GSM922201_Sample_G9.expression.txt",
  "GSM896810_Sample_B2.expression.txt", "GSM922202_Sample_G10.expression.txt",
  "GSM896811_Sample_B3.expression.txt", "GSM922203_Sample_G11.expression.txt",
  "GSM896812_Sample_B4.expression.txt", "GSM922204_Sample_G12.expression.txt",
  "GSM896813_Sample_B5.expression.txt", "GSM922205_Sample_G13.expression.txt",
  "GSM896814_Sample_B6.expression.txt", "GSM922206_Sample_H1.expression.txt",
  "GSM896815_Sample_M1.expression.txt", "GSM922207_Sample_H2.expression.txt",
  "GSM896816_Sample_M2.expression.txt", "GSM922208_Sample_H3.expression.txt",
  "GSM922146_Sample_C4.expression.txt", "GSM922209_Sample_H4.expression.txt",
  "GSM922147_Sample_C5.expression.txt", "GSM922210_Sample_H5.expression.txt",
  "GSM922148_Sample_C6.expression.txt", "GSM922211_Sample_H6.expression.txt",
  "GSM922149_Sample_C7.expression.txt", "GSM922212_Sample_H7.expression.txt",
  "GSM922150_Sample_L1.expression.txt", "GSM922213_Sample_H8.expression.txt",
  "GSM922151_Sample_L2.expression.txt", "GSM922214_Sample_H9.expression.txt",
  "GSM922152_Sample_L3.expression.txt", "GSM922215_Sample_H10.expression.txt",
  "GSM922153_Sample_L4.expression.txt", "GSM922216_Sample_R1.expression.txt",
  "GSM922154_Sample_L5.expression.txt", "GSM922217_Sample_R2.expression.txt",
  "GSM922155_Sample_L6.expression.txt", "GSM922218_Sample_R3.expression.txt",
  "GSM922156_Sample_L7.expression.txt", "GSM922219_Sample_R4.expression.txt",
  "GSM922157_Sample_L8.expression.txt", "GSM922220_Sample_R5.expression.txt",
  "GSM922158_Sample_A1.expression.txt", "GSM922221_Sample_R6.expression.txt",
  "GSM922159_Sample_A2.expression.txt", "GSM922222_Sample_R7.expression.txt",
  "GSM922160_Sample_A3.expression.txt", "GSM922223_Sample_R8.expression.txt",
  "GSM922161_Sample_A4.expression.txt", "GSM922224_Sample_M4.expression.txt",
  "GSM922162_Sample_A5.expression.txt", "GSM922225_Sample_M5.expression.txt",
  "GSM922163_Sample_A6.expression.txt", "GSM922226_Sample_M6.expression.txt",
  "GSM922164_Sample_A7.expression.txt", "GSM922227_Sample_M7.expression.txt",
  "GSM922165_Sample_A8.expression.txt", "GSM922228_Sample_M8.expression.txt",
  "GSM922166_Sample_E1.expression.txt", "GSM922230_Sample_M10.expression.txt",
  "GSM922167_Sample_E2.expression.txt", "GSM922250_Sample_P1.expression.txt",
  "GSM922168_Sample_E3.expression.txt", "GSM922251_Sample_P2.expression.txt",
  "GSM922169_Sample_E4.expression.txt", "GSM922252_Sample_P5.expression.txt",
  "GSM922170_Sample_K1.expression.txt", "GSM922253_Sample_P7.expression.txt",
  "GSM922171_Sample_K2.expression.txt", "GSM922254_Sample_P9.expression.txt",
  "GSM922172_Sample_K3.expression.txt", "GSM922255_Sample_P11.expression.txt",
  "GSM922173_Sample_K4.expression.txt", "GSM922256_Sample_P12.expression.txt",
  "GSM922174_Sample_K5.expression.txt", "GSM922257_Sample_P13.expression.txt",
  "GSM922175_Sample_K6.expression.txt", "GSM922258_Sample_P14.expression.txt",
  "GSM922176_Sample_K7.expression.txt", "GSM922259_Sample_P15.expression.txt",
  "GSM922177_Sample_K8.expression.txt", "GSM922260_Sample_P17.expression.txt",
  "GSM922178_Sample_D1.expression.txt", "GSM922261_Sample_S1.expression.txt",
  "GSM922179_Sample_D2.expression.txt", "GSM922262_Sample_S2.expression.txt",
  "GSM922180_Sample_D3.expression.txt", "GSM922263_Sample_S3.expression.txt",
  "GSM922181_Sample_D4.expression.txt", "GSM922264_Sample_S4.expression.txt",
  "GSM922182_Sample_D5.expression.txt", "GSM922265_Sample_S5.expression.txt",
  "GSM922183_Sample_D6.expression.txt", "GSM922266_Sample_S6.expression.txt",
  "GSM922184_Sample_D7.expression.txt", "GSM922267_Sample_S7.expression.txt",
  "GSM922185_Sample_D8.expression.txt", "GSM922268_Sample_S8.expression.txt",
  "GSM922186_Sample_F1.expression.txt", "GSM922269_Sample_S9.expression.txt",
  "GSM922187_Sample_F2.expression.txt", "GSM922270_Sample_S10.expression.txt",
  "GSM922188_Sample_F3.expression.txt", "GSM922271_Sample_S11.expression.txt",
  "GSM922189_Sample_F4.expression.txt", "GSM922272_Sample_S12.expression.txt",
  "GSM922190_Sample_F5.expression.txt", "GSM922273_Sample_S13.expression.txt",
  "GSM922191_Sample_F6.expression.txt", "GSM922274_Sample_S14.expression.txt",
  "GSM922192_Sample_F7.expression.txt", "GSM922275_Sample_S15.expression.txt",
  "GSM922193_Sample_F8.expression.txt", "GSM922194_Sample_G2.expression.txt"
)

# Initialize an empty list to store individual data frames
count_data_list <- list()

# Read each file and store its data frame in the list
for (filepath in file_paths) {
  count_data <- read_expression_file(filepath)
  sample_name <- sub(".*_(Sample_.*).expression.txt", "\\1", filepath)
  colnames(count_data)[2] <- sample_name
  count_data_list[[sample_name]] <- count_data
}

# Combine all data frames by the 'gene_name' column
combined_counts <- Reduce(function(x, y) merge(x, y, by = "gene_name", all = TRUE), count_data_list)

# Replace NAs with 0
combined_counts[is.na(combined_counts)] <- 0
