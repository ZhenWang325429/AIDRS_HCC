#########################################
######### Data Integration Section ##########
#########################################

# Load necessary libraries
library(tidyverse)           # For data manipulation
library(here)                # For handling file paths
library(SummarizedExperiment) # For managing experiment data

# Load previously processed data
load(here('D:/2023/MOVICS_HCC/Data/', 'Processed_mRNA_Data.Rdata'))
load(here('D:/2023/MOVICS_HCC/Data/LIHC_Rawdata', 'TCGA-LIHC_CNV.Rdata'))
load(here('D:/2023/MOVICS_HCC/Data/LIHC_Rawdata', 'TCGA-LIHC_METHY_beta.Rdata'))
load(here('D:/2023/MOVICS_HCC/Data/LIHC_Rawdata', 'TCGA-LIHC_SNP.Rdata'))
load(here('D:/2023/MOVICS_HCC/Data/LIHC_Rawdata', 'Survival_SupplementalTable_S1_20171025_xena_sp'))

# Clinical data processing
Pd <- vroom(here('D:/2023/MOVICS_HCC/Data/LIHC_Rawdata', 'Survival_SupplementalTable_S1_20171025_xena_sp')) %>%
  filter(`cancer type abbreviation` == 'LIHC') %>%
  filter(sample %in% (colnames(mRNA_expr_counts_symbol_max))) %>%
  mutate(ID = sample) %>%
  select("ID", "OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time") %>%
  mutate(bcr_patient_barcode = substr(ID, 1, 12))

# Intersect all omics data and clinical data
intersect_result <- Reduce(intersect, list(colnames(mRNA_expr_counts_symbol_max),
                                           colnames(Meth_beta),
                                           colnames(SNP_matrix_0_1),
                                           unique(CNV_segment$sample),
                                           unique(SNP_maf$Tumor_Sample_Barcode),
                                           Pd_Cli$ID))

# Subset data by intersected sample IDs
mRNA_expr_counts_symbol_max <- mRNA_expr_counts_symbol_max %>%
  select(intersect_result)
mRNA_expr_fpkm_symbol_max <- mRNA_expr_fpkm_symbol_max %>%
  select(intersect_result)
mRNA_expr_tpm_symbol_max <- mRNA_expr_tpm_symbol_max %>%
  select(intersect_result)

# Similarly, subset other data (miRNA, lncRNA, Meth, SNP) by intersected sample IDs
# Save integrated data
save(mRNA_expr_counts_symbol_max, mRNA_expr_fpkm_symbol_max, mRNA_expr_tpm_symbol_max,
     Meth_beta, SNP_matrix_0_1, Pd_Cli,
     file = here('D:/2023/MOVICS_HCC/Data/', 'Integrated_Data.Rdata'))

# Run survival analysis using Cox regression, select relevant omics features
surv.info <- Pd_Cli

# Example: Select top features from mRNA data
elite.mRNA <- getElites(dat = log2(mRNA_expr_tpm_symbol_max + 1),
                        method = "cox",
                        surv.info = surv.info,
                        p.cutoff = 0.001,
                        elite.num = 500)

# Filter top features
top_500 <- elite.mRNA$unicox.res %>%
  arrange(pvalue) %>%
  select(gene) %>%
  head(500)

mRNA <- elite.mRNA$elite.dat %>%
  as.data.frame() %>%
  filter(rownames(.) %in% top_500$gene)

# Save elite features for further analysis
save(elite.mRNA, elite.miRNA, elite.lncRNA,
     elite.SNP, elite.Meth, CNV_segment, SNP_maf,
     file = here('D:/2023/MOVICS_HCC/Data/', 'Elite_Features.Rdata'))

# Perform multi-omics integrative clustering and save results
optk.LIHC <- getClustNum(data = LIHC_data[1:5],
                         is.binary = c(F, F, F, F, T),
                         try.N.clust = 2:8,
                         fig.name = "CLUSTER NUMBER OF TCGA-LIHC",
                         fig.path = 'D:/2023/MOVICS_HCC/Results/')

# Clustering and consensus analysis
moic.res.list <- getMOIC(data = LIHC_data[1:5],
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"),
                         N.clust = 3,
                         type = c("gaussian", "gaussian", "gaussian", "gaussian", "binomial"))

# Save clustering results
save(iClusterBayes.res, moic.res.list, cmoic.LIHC,
     file = here('D:/2023/MOVICS_HCC/Data/', 'Cluster_Results.Rdata'))
