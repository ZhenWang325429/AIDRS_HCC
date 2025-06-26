
####Drugs_Sensibility#####


# Load necessary libraries
library(limma)
library(oncoPredict)
library(parallel)
library(ggplot2)
library(ggpubr)
library(future.apply)
library(car)
library(ridge)
library(rstatix)
library(vroom)
library(dplyr)
library(tidyr)

# Set seed for reproducibility
set.seed(12345)

# Define p-value filter threshold
pFilter = 0.05

# Load the required datasets (paths can be customized)
data_dir <- 'F:/2023/'
load(paste0(data_dir, 'Disulfidptosis/2_Filter_Data.Rdata'))
load(paste0(data_dir, 'MOVICS_HCC/MOVIC_HCC/Data_V4/06_NTP_Other_Sets.Rdata'))
load(paste0(data_dir, 'MOVICS_HCC/MOVIC_HCC/Data_V4/03_LIHC_data_Muti_Omics.Rdata'))

# **Multi-threaded drug sensitivity analysis**

# Convert mRNA expression matrix to numeric
dat <- apply(mRNA_expr_counts_symbol_max, 2, as.numeric)

# Set rownames for expression matrix
rownames(dat) <- rownames(mRNA_expr_counts_symbol_max)

# Create lists for data matrices and predictions
dat_list <- list(
  GSE141200 = Exp_GSE141200_Counts_filter,
  GSE14520 = Exp_GSE14520_Array_filter,
  ICGC = Exp_ICGC_Counts_filter,
  TCGA = dat
)

ntp_list <- list(
  GSE141200 = NTP.pred_GSE141200,
  GSE14520 = NTP.pred_GSE14520,
  ICGC = NTP.pred_ICGC,
  TCGA = cmoic.LIHC
)

# Drug data and statistics lists
Drug_dat <- list()
Drug_dat_Sta <- list()

# Iterate over the datasets
for (i in 1:length(dat_list)) {
  
  # Read and process drug sensitivity predictions
  Drug_dat[[i]] <- vroom::vroom(file = paste0(data_dir, 'MOVICS_HCC/MOVIC_HCC/Results_V4/calcPhenotype_Output/', IDs[i], "_DrugPredictions.csv")) %>% 
    column_to_rownames("...1") %>% 
    select(colnames(.)[!(duplicated(str_remove(colnames(.),"_\\d+")))]) %>% 
    rename_at(vars(matches("_\\d+")), ~str_remove(., "_\\d+")) %>% 
    mutate(Group = ifelse(ntp_list[[i]]$clust.res$clust == 1, 'CS1', 'CS2') %>% factor(., levels = c('CS1', 'CS2'))) %>% 
    rownames_to_column('ID') %>% 
    pivot_longer(cols = -c('ID','Group'), names_to = 'Drugs', values_to = 'Score')
  
  # Calculate drug sensitivity statistics and adjust p-values
  Drug_dat_Sta[[i]] <- Drug_dat[[i]] %>%
    group_by(Drugs) %>%
    na.omit() %>% 
    filter(!Score == 'Inf') %>% 
    mutate(Score = log2(Score + 1)) %>% 
    wilcox_test(Score ~ Group, comparisons = c('CS1', 'CS2'), detailed = TRUE) %>%
    adjust_pvalue(method = 'BH') %>%
    add_significance("p.adj") %>%
    arrange(p)
}

# Identify drugs with significant changes across datasets
Drug_Up_ids <- intersect(Drug_dat_Sta[[3]][Drug_dat_Sta[[3]]$estimate > 0 & Drug_dat_Sta[[3]]$p < 0.05,]$Drugs,
                         intersect(Drug_dat_Sta[[2]][Drug_dat_Sta[[2]]$estimate > 0 & Drug_dat_Sta[[2]]$p < 0.05,]$Drugs,
                                   Drug_dat_Sta[[4]][Drug_dat_Sta[[4]]$estimate > 0 & Drug_dat_Sta[[4]]$p < 0.05,]$Drugs))

Drug_Down_ids <- intersect(Drug_dat_Sta[[3]][Drug_dat_Sta[[3]]$estimate < 0 & Drug_dat_Sta[[3]]$p < 0.05,]$Drugs,
                           intersect(Drug_dat_Sta[[2]][Drug_dat_Sta[[2]]$estimate < 0 & Drug_dat_Sta[[2]]$p < 0.05,]$Drugs,
                                     Drug_dat_Sta[[4]][Drug_dat_Sta[[4]]$estimate < 0 & Drug_dat_Sta[[4]]$p < 0.05,]$Drugs))

# **Plot heatmap to visualize results**

# Prepare data for heatmap
dat1 <- Drug_dat_Sta[[3]] %>% 
  as.data.frame() %>% 
  filter(Drugs %in% c(Drug_Up_ids, Drug_Down_ids)) %>% 
  arrange(estimate) %>% 
  mutate(Drugs = factor(c(Drug_Up_ids, Drug_Down_ids), levels = c(Drug_Up_ids, Drug_Down_ids))) %>% 
  mutate(ID = 'ICGC-LIRI') %>% 
  mutate(Label = ifelse(p < 0.001, '***', ifelse(p < 0.01, '**', '*')))

dat2 <- Drug_dat_Sta[[2]] %>% 
  as.data.frame() %>% 
  filter(Drugs %in% c(Drug_Up_ids, Drug_Down_ids)) %>%  
  arrange(estimate) %>% 
  mutate(Drugs = factor(c(Drug_Up_ids, Drug_Down_ids), levels = c(Drug_Up_ids, Drug_Down_ids))) %>% 
  mutate(ID = 'GSE14520') %>% 
  mutate(Label = ifelse(p < 0.001, '***', ifelse(p < 0.01, '**', '*')))

dat3 <- Drug_dat_Sta[[4]] %>%
  as.data.frame() %>% 
  filter(Drugs %in% c(Drug_Up_ids, Drug_Down_ids)) %>%  
  arrange(estimate) %>% 
  mutate(Drugs = factor(c(Drug_Up_ids, Drug_Down_ids), levels = c(Drug_Up_ids, Drug_Down_ids))) %>% 
  mutate(ID = 'TCGA-LIHC') %>% 
  mutate(Label = ifelse(p < 0.001, '***', ifelse(p < 0.01, '**', '*')))

# Function to scale data to the range [-1, 1]
scale_to_minus1_1 <- function(x) {
  x_non_na <- na.omit(x)
  if (length(x_non_na) == 0 || var(x_non_na) == 0) {
    return(rep(NA, length(x)))  # If data is constant or empty, return NA
  } else {
    x_scaled <- (x - mean(x_non_na)) / sd(x_non_na)
    return(2 * (x_scaled - min(x_scaled)) / (max(x_scaled) - min(x_scaled)) - 1)
  }
}

# Prepare heatmap data
Heat_dat <- data.frame('Drugs' = dat1$Drugs,
                       dat3 = dat3$estimate,
                       dat1 = dat1$estimate,
                       dat2 = dat2$estimate,
                       dat3_label = dat3$Label,
                       dat1_label = dat1$Label,
                       dat2_label = dat2$Label) %>% 
  mutate('TCGA-LIHC' = scale_to_minus1_1(dat3),
         'ICGC-LIRI' = scale_to_minus1_1(dat1),
         'GSE14520' = scale_to_minus1_1(dat2)) %>% 
  column_to_rownames('Drugs')

# Reshape data for heatmap plotting
heat_dat <- Heat_dat %>% 
  rownames_to_column("Drugs") %>% 
  select(c("Drugs", "TCGA-LIHC", "ICGC-LIRI", "GSE14520", "dat3_label", "dat1_label", "dat2_label")) %>%
  pivot_longer(cols = c('TCGA-LIHC', 'ICGC-LIRI', 'GSE14520'), names_to = "dataset", values_to = "score") %>%
  pivot_longer(cols = c('dat3_label', 'dat1_label', 'dat2_label'), names_to = "label_set", values_to = "label") %>%
  filter((dataset == "TCGA-LIHC" & label_set == "dat3_label") |
           (dataset == "ICGC-LIRI" & label_set == "dat1_label") |
           (dataset == "GSE14520" & label_set == "dat2_label")) %>%
  select(-label_set)

# Reorder Drugs factor
heat_dat$Drugs <- factor(heat_dat$Drugs, levels = c(Drug_Up_ids, Drug_Down_ids))

# Plot heatmap using ggplot2
Drugs_heatmap <- ggplot(heat_dat, aes(x = dataset, y = Drugs, fill = score)) + 
  scale_x_discrete(position = 'bottom') + 
  scale_fill_gradient2(low = "#4AA4BE", high = "#F0A13C", mid = "white", midpoint = 0.25) +
  geom_tile(width = 1, height = 1) +
  theme_minimal() +
  geom_text(aes(label = label), color = "black", size = 8) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 16),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 20, face = 'bold'),
        legend.position = 'bottom',
        panel.grid.major = element_blank()) +
  geom_hline(yintercept = seq(from = 1.5, to = 66.5, by = 1), size = .8, color = 'white') +
  geom_vline(xintercept = c(1.5, 2.5), size = .8, color = 'white') +
  xlab(NULL) + ylab(NULL)

# Preview and save heatmap plot
ggpreview(Drugs_heatmap, width = 5, height = 20)
ggsave(filename = paste0(data_dir, 'MOVICS_HCC/MOVIC_HCC/Results_V4/Drugs_heatmap.pdf'), width = 5, height = 20)
