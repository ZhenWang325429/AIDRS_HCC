######################################
######### Data Download Section ##########
######################################

# Required Libraries
library(MOVICS)             # For multi-omics integration and analysis
library(tidyverse)          # For data wrangling
library(vroom)              # For fast data loading
library(here)               # For creating file paths
library(SummarizedExperiment) # For managing experiment data

# Load raw data
# Customize the file paths as per your directory structure
mRNA_data_path <- here('D:/2023/MOVICS_HCC/Data/LIHC_Rawdata', 'TCGA-LIHC_mRNA.Rdata')
load(mRNA_data_path)

# Extract RNA data: Counts, TPM, and FPKM matrices
RNA_expr_counts <- assay(data, "unstranded")
RNA_expr_tpm <- assay(data, "tpm_unstrand")
RNA_expr_fpkm <- assay(data, "fpkm_unstrand")

# Extract gene symbols
symbol <- rowData(data)$gene_name

# Bind gene symbols with expression data
RNA_expr_counts_symbol <- cbind(data.frame(symbol), as.data.frame(RNA_expr_counts))
RNA_expr_tpm_symbol <- cbind(data.frame(symbol), as.data.frame(RNA_expr_tpm))
RNA_expr_fpkm_symbol <- cbind(data.frame(symbol), as.data.frame(RNA_expr_fpkm))

# Filter for protein-coding genes and remove duplicates
mRNA_expr_counts_symbol_max <- RNA_expr_counts_symbol %>%
  as_tibble() %>%
  filter(data@rowRanges$gene_type == 'protein_coding') %>%
  mutate(meanrow = rowMeans(.[,-1]), .before = 2) %>%
  arrange(desc(meanrow)) %>%
  distinct(symbol, .keep_all = T) %>%
  select(-meanrow) %>%
  column_to_rownames(var = "symbol") %>%
  select(-which(duplicated(substr(names(.), 1, 15)))) %>%
  rename_all(~substr(., 1, 15)) %>%
  select(-which(substr(names(.), 14, 15) == "11"))

# Repeat for TPM and FPKM data
mRNA_expr_tpm_symbol_max <- RNA_expr_tpm_symbol %>%
  as_tibble() %>%
  filter(data@rowRanges$gene_type == 'protein_coding') %>%
  mutate(meanrow = rowMeans(.[,-1]), .before = 2) %>%
  arrange(desc(meanrow)) %>%
  distinct(symbol, .keep_all = T) %>%
  select(-meanrow) %>%
  column_to_rownames(var = "symbol") %>%
  select(-which(duplicated(substr(names(.), 1, 15)))) %>%
  rename_all(~substr(., 1, 15)) %>%
  select(-which(substr(names(.), 14, 15) == "11"))

mRNA_expr_fpkm_symbol_max <- RNA_expr_fpkm_symbol %>%
  as_tibble() %>%
  filter(data@rowRanges$gene_type == 'protein_coding') %>%
  mutate(meanrow = rowMeans(.[,-1]), .before = 2) %>%
  arrange(desc(meanrow)) %>%
  distinct(symbol, .keep_all = T) %>%
  select(-meanrow) %>%
  column_to_rownames(var = "symbol") %>%
  select(-which(duplicated(substr(names(.), 1, 15)))) %>%
  rename_all(~substr(., 1, 15)) %>%
  select(-which(substr(names(.), 14, 15) == "11"))

# Save the processed mRNA expression matrices
save(mRNA_expr_counts_symbol_max, mRNA_expr_fpkm_symbol_max, mRNA_expr_tpm_symbol_max,
     file = here('D:/2023/MOVICS_HCC/Data/', 'Processed_mRNA_Data.Rdata'))

# Similarly, repeat for miRNA, lncRNA, CNV, Meth, SNP as in the original code...

