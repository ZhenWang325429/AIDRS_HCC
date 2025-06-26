###TIDE#######

# Load required libraries
library(vroom)
library(tidyverse)

# Step 1: Load the TCGA-LIHC data
# Load expression data
load('F:/2023/MOVICS_HCC/MOVIC_HCC/Data/LIHC_Rawdata/TCGA-LIHC_mRNA.Rdata')
load('F:/2023/MOVICS_HCC/MOVIC_HCC/Data_V4/03_LIHC_data_Muti_Omics.Rdata')

# Prepare Counts matrix for analysis
LIHC_Counts <- data@assays@data$unstranded  # Non-log transformed counts
dim(LIHC_Counts)

# Add column and row names
colnames(LIHC_Counts) <- str_sub(data@colData@rownames, 1, 15)
rownames(LIHC_Counts) <- data@rowRanges$gene_name

# Filter for cancer samples
LIHC_Counts <- LIHC_Counts[, as.numeric(str_sub(colnames(LIHC_Counts), 14, 15)) < 10]
LIHC_Counts <- LIHC_Counts[, str_sub(colnames(LIHC_Counts), 1, 15) %in% colnames(LIHC_data$mRNA)]
LIHC_Counts <- LIHC_Counts[, !duplicated(colnames(LIHC_Counts))]

# Prepare TPM matrix
LIHC_TPM <- data@assays@data$tpm_unstrand  # Non-log transformed TPM
dim(LIHC_TPM)

# Add column and row names
colnames(LIHC_TPM) <- str_sub(data@colData@rownames, 1, 15)
rownames(LIHC_TPM) <- data@rowRanges$gene_name

# Filter for samples present in the counts matrix
LIHC_TPM <- LIHC_TPM[, colnames(LIHC_TPM) %in% colnames(LIHC_Counts)]

# Standardize TPM matrix
mydat <- t(apply(LIHC_TPM, 1, function(x) x - mean(x)))

# Save the standardized TPM matrix for further processing
write.table(mydat, file = './Data_V4/TCGA-LIHC_TIDE.txt', sep = '\t')

# Load TIDE results (predicted immune response scores)
LIHC_TIDE_Res <- vroom('F:/2023/MOVICS_HCC/Data_V4/TIDE/TCGA_LIHC_res.txt') %>%
  column_to_rownames('...1') %>%
  filter(rownames(.) %in% cmoic.LIHC$clust.res$samID)

# Check the responder distribution
table(LIHC_TIDE_Res$Responder)  # Count TRUE and FALSE responders

# Reorder TIDE results to match the sample IDs from the clustering results
LIHC_TIDE_Res = LIHC_TIDE_Res[match(cmoic.LIHC$clust.res$samID, rownames(LIHC_TIDE_Res)),]

# Save the TIDE results
save(LIHC_TIDE_Res, file = 'F:/2023/MOVICS_HCC/Data_V4/TIDE/TCGA_LIHC_TIDE_Res.Rdata')
