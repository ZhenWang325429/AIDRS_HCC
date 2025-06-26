# SNV Analysis ##########

# Load required libraries
library(TCGAbiolinks)
library(maftools)

# (1) TMB and MATH Heterogeneity Calculation -------------------------------------

# Load MAF data
load('F:/2023/MOVICS_HCC/MOVIC_HCC/Data/LIHC_Rawdata/TCGA-LIHC_SNP.Rdata')

# Remove "chr" prefix and reformat chromosome
maf@data <- maf@data %>%
  mutate(Chromosome = str_remove(Chromosome, "chr")) %>%
  mutate(Chromosome = ifelse(Chromosome == 'X', 23, Chromosome)) %>%
  mutate(Chromosome = as.numeric(Chromosome))

# Calculate TMB and MATH
TMB <- tmb(maf = maf, captureSize = 50, logScale = TRUE)

# Visualize TMB differences between CS1 and CS2
ggsave("F:/2023/MOVICS_HCC/MOVIC_HCC/Results_V4/TME/TMB.pdf", p, width = 4, height = 4)

# (2) Mutation Landscape and Oncoplot ---------------------------------------------

# Create oncoplot for mutation data
maf@clinical.data <- maf@clinical.data %>%
  mutate(Cluster = ifelse(str_sub(Tumor_Sample_Barcode, 1, 15) %in% cmoic.LIHC$clust.res[cmoic.LIHC$clust.res$clust == 1,]$samID, "CS1", "CS2"))

pdf(file = 'F:/2023/MOVICS_HCC/MOVIC_HCC/Results_V4/Oncoplot.pdf', width = 12, height = 6)
oncoplot(maf = maf, top = 20, drawColBar = FALSE, clinicalFeatures = "Cluster")
dev.off()

# (3) Mutation Frequency Comparison ---------------------------------------------

# Compare mutation frequencies for TP53 and CTNNB1 in CS1 and CS2
table(unique(str_sub(data[data$Hugo_Symbol == 'TP53',]$Tumor_Sample_Barcode, 1, 15)) %in% cmoic.LIHC$clust.res[cmoic.LIHC$clust.res$clust == 1,]$samID)

# Plot mutation proportions for TP53 and CTNNB1
ggsave(filename = 'F:/2023/MOVICS_HCC/MOVIC_HCC/Results_V4/TP53_mut.pdf', p_TP53, width = 3, height = 4)
ggsave(filename = 'F:/2023/MOVICS_HCC/MOVIC_HCC/Results_V4/CTNNB1_mut.pdf', p_CTNNB1, width = 3, height = 4)

# (4) Identify Significant Mutations Between Groups -----------------------------

# Perform enrichment analysis for mutations in CS1 and CS2
fab.ce <- clinicalEnrichment(maf = maf, clinicalFeature = 'Cluster')

# Plot enrichment results
plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05)

# Save results
save(fab.ce, file = "F:/2023/MOVICS_HCC/MOVIC_HCC/Results_V4/TME/TCGA_LIHC_SNV_Res.Rdata")
