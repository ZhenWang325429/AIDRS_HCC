# CNV_Analysis#########

# Load required libraries
library(TCGAbiolinks)
library(tidyverse)
library(vroom)

# (1) Load data and prepare segment and marker files --------------------------------

# Load raw CNV data and clustering results
load('F:/2023/MOVICS_HCC/MOVIC_HCC/Data/LIHC_Rawdata/TCGA-LIHC_CNV.Rdata')
load('F:/2023/MOVICS_HCC/MOVIC_HCC/Data_V4/04_LIHC_moic.res.list.Rdata')

# Inspect the first few rows of CNV data
head(data)

# Filter tumor (1~9) and non-tumor (>=10) samples
table(str_sub(data$Sample,14,15) >= 10)

# Prepare CNV data for CS1 and CS2 clusters
LIHC_CNV_data <- data %>%
  filter(str_sub(Sample,14,15) < 10) %>%
  select("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean") %>%
  rename(`Start Position` = Start, `End Position` = End, `Num markers` = Num_Probes, Seg.CN = Segment_Mean) %>%
  mutate(Cluster = ifelse(str_sub(Sample,1,15) %in% rownames(cmoic.LIHC$clust.res[cmoic.LIHC$clust.res$clust == 1,]), "CS1", "CS2"))

# Output filtered data for CS1 and CS2 clusters
write.table(LIHC_CNV_data %>% filter(Cluster == "CS1"), file = "F:/2023/MOVICS_HCC/MOVIC_HCC/Data_V4/TCGA-LIHC_GISTIC2_Input/LIHC_segment_file_CS1.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(LIHC_CNV_data %>% filter(Cluster == "CS2"), file = "F:/2023/MOVICS_HCC/MOVIC_HCC/Data_V4/TCGA-LIHC_GISTIC2_Input/LIHC_segment_file_CS2.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

# (2) Prepare marker file ------------------------------------------------------------

# Load marker file and filter
hg38_marker_file <- read.delim("F:/2023/MOVICS_HCC/MOVIC_HCC/Data_V4/TCGA-LIHC_GISTIC2_Input/snp6.na35.remap.hg38.subset.txt.gz")
hg_marker_file <- hg38_marker_file %>%
  filter(freqcnv == "FALSE") %>%
  select(1:3) %>%
  rename(`Marker Name` = probeid, Chromosome = chr, `Marker Position` = pos)

# Save marker file
write.table(hg_marker_file, file = "F:/2023/MOVICS_HCC/MOVIC_HCC/Data_V4/TCGA-LIHC_GISTIC2_Input/hg38_marker_file.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

# (3) CNV Distribution: Histogram and Donut Plot -----------------------------------

# Load CNV scores for CS1 and CS2
scores_CS1 <- read.table("F:/2023/MOVICS_HCC/MOVIC_HCC/Data_V4/TCGA-LIHC_GISTIC2_Res/CS1_Res/scores.gistic", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  as.data.frame() %>%
  mutate(Cluster = "CS1")
scores_CS2 <- read.table("F:/2023/MOVICS_HCC/MOVIC_HCC/Data_V4/TCGA-LIHC_GISTIC2_Res/CS2_Res/scores.gistic", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  as.data.frame() %>%
  mutate(Cluster = "CS2")

# Combine CS1 and CS2 data
scores <- rbind(scores_CS2, scores_CS1)

# Plot CNV distribution as histogram and donut plot
library(ggplot2)
p <- ggplot(scores, aes(StartPos, Frequency)) +
  geom_area(aes(group = Type, fill = Type)) +
  scale_fill_manual(values = c("#a9abcb", "#d49c9c")) +
  theme_classic() +
  facet_grid(Cluster ~ .)

# Save CNV distribution plot
ggsave("F:/2023/MOVICS_HCC/MOVIC_HCC/Results_V4/TME/CNV_Dif_Frequency.pdf", plot = p, width = 14, height = 4)

# (4) Further CNV-Gene Expression Analysis ----------------------------------------

# Filter and prepare CNV data for CS1 and CS2
CS1_CNV <- read.table("F:/2023/MOVICS_HCC/MOVIC_HCC/Data_V4/TCGA-LIHC_GISTIC2_Res/CS1_Res/broad_data_by_genes.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  as.data.frame() %>%
  column_to_rownames('Gene.Symbol') %>%
  select(-c(1:2)) %>%
  rename_with(~ str_replace_all(., "\\.", "-")) %>%
  pivot_longer(cols = 2:ncol(.), values_to = 'CNV', names_to = 'ID') %>%
  mutate(Cluster = "CS1")

CS2_CNV <- read.table("F:/2023/MOVICS_HCC/MOVIC_HCC/Data_V4/TCGA-LIHC_GISTIC2_Res/CS2_Res/broad_data_by_genes.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  as.data.frame() %>%
  column_to_rownames('Gene.Symbol') %>%
  select(-c(1:2)) %>%
  rename_with(~ str_replace_all(., "\\.", "-")) %>%
  pivot_longer(cols = 2:ncol(.), values_to = 'CNV', names_to = 'ID') %>%
  mutate(Cluster = "CS2")

# Perform differential analysis using DESeq2 for CNV and gene expression
library(DESeq2)
load('F:/2023/MOVICS_HCC/MOVIC_HCC/Data_V4/03_LIHC_data_Muti_Omics.Rdata')

# Differential analysis of CNV and gene expression
res <- DESeq2::results(dds)

# Plot boxplots for gene expression and CNV comparisons between CS1 and CS2
ggsave("F:/2023/MOVICS_HCC/MOVIC_HCC/Results_V4/TME/TPM_CS.pdf", plot = p1, width = 4, height = 3)
ggsave("F:/2023/MOVICS_HCC/MOVIC_HCC/Results_V4/TME/CNV_CS.pdf", plot = p2, width = 4, height = 3)

# Save CNV and gene expression data results
save(CNV_TPM, file = "F:/2023/MOVICS_HCC/MOVIC_HCC/Results_V4/TME/TCGA_LIHC_CNV_Res.Rdata")
