######## Clinical_Feature #################

# Load the required data
load(here('F:/2023/MOVICS_HCC/MOVIC_HCC/Data/','01_Data_Muti_Omics.Rdata'))
load(here('F:/2023/MOVICS_HCC/MOVIC_HCC/Data_V4/','04_LIHC_moic.res.list.Rdata'))
load(here('F:/2023/MOVICS_HCC/MOVIC_HCC/','Data_V4/06_NTP_Other_Sets.Rdata'))

# **Prepare TCGA Clinical Feature Data**
Cli_data <- Pd_Cli %>%
  mutate(Cluster = ifelse(ID %in% cmoic.LIHC$clust.res[(cmoic.LIHC$clust.res)$clust == 1,]$samID, "CS1", "CS2") %>%
           factor(levels = c("CS1", "CS2")),
         Age = ifelse(Age == ">60", ">=60", "<60") %>% factor(levels = c("<60", ">=60")),
         Race = ifelse(Race == "Asian", "Asian", ifelse(Race == "White", "White", NA)) %>% factor(levels = c('Asian', 'White')),
         Child = ifelse(Child %in% c('B', 'C'), 'B/C', 'A') %>% factor(levels = c('A', 'B/C')),
         Grade = ifelse(Grade %in% c('G1', 'G2'), 'G1/G2', 'G3/G4') %>% factor(levels = c('G1/G2', 'G3/G4')),
         Stage = ifelse(Stage %in% c('Stage I', 'Stage II'), 'I/II', 'III/IV') %>% factor(levels = c('I/II', 'III/IV')),
         PT = ifelse(PT %in% c('<10', '10-14'), '<14', '>=14') %>% factor(levels = c('<14', '>=14')),
         ALB = ifelse(ALB %in% c('<3.5', '3.5-5.5'), '<5.5', '>=5.5') %>% factor(levels = c('<5.5', '>=5.5')),
         Time = OS.time, Status = OS) %>%
  select(c("ID", "Time", "Status", "Cluster", "Age", "Gender", "Race", "Event", "Child", "Grade", "Stage", "Hepatitis", "Alcohol", "AFP", "PT", "ALB")) %>%
  arrange(Cluster)

# **Statistical analysis of clinical features and Cluster**
com_res <- Cli_data %>% summarise(
  p_Age = chisq.test(table(Age, Cluster))$p.value,
  p_Gender = chisq.test(table(Gender, Cluster))$p.value,
  p_Race = chisq.test(table(Race, Cluster))$p.value,
  p_Event = chisq.test(table(Event, Cluster))$p.value,
  p_Child = chisq.test(table(Child, Cluster))$p.value,
  p_Grade = chisq.test(table(Grade, Cluster))$p.value,
  p_Stage = chisq.test(table(Stage, Cluster))$p.value,
  p_Hepatitis = chisq.test(table(Hepatitis, Cluster))$p.value,
  p_Alcohol = chisq.test(table(Alcohol, Cluster))$p.value,
  p_AFP = chisq.test(table(AFP, Cluster))$p.value,
  p_PT = chisq.test(table(PT, Cluster))$p.value,
  p_ALB = chisq.test(table(ALB, Cluster))$p.value
)

# **Plot heatmap with clinical features annotations**
library(ComplexHeatmap)
ha <- HeatmapAnnotation(df = Cli_data[, 4:16],
                        col = list(Cluster = c('CS1' = "#2EC4B6", 'CS2' = "#E71D36"),
                                   Age = c('<60' = '#4DAF4A', '>=60' = '#984EA3'),
                                   Gender = c('Male' = '#FFFFBF', 'Female' = '#A65628'),
                                   Race = c('Asian' = '#BC80BD', 'White' = '#BEBAD5'),
                                   Event = c('Alive' = '#FC7300', 'Dead' = '#F781BF'),
                                   Child = c('A' = '#BF5B17', 'B/C' = '#ABDDA4'),
                                   Grade = c('G1/G2' = '#FDC086', 'G3/G4' = '#BEAED4'),
                                   Stage = c('I/II' = '#E6F598', 'III/IV' = '#FDAE61'),
                                   Hepatitis = c('Yes' = 'black', 'No' = '#898B79'),
                                   Alcohol = c('Yes' = 'black', 'No' = '#898B79'),
                                   AFP = c('<200' = '#8DD3C7', '>=200' = '#B3DE69'),
                                   PT = c('<14' = '#5A8AA2', '>=14' = '#AC726D'),
                                   ALB = c('<5.5' = '#F090A2', '>=5.5' = '#466780')),
                        annotation_name_side = 'left')

# Create heatmap for clinical features
zero_row_mat <- matrix(nrow = 0, ncol = nrow(Cli_data))
Hm <- Heatmap(zero_row_mat,
              row_title_side = 'left',
              top_annotation = ha)

# Save heatmap as PDF
pdf(file = "F:/2023/MOVICS_HCC/Results/TCGA_Cli.pdf", width = 10, height = 6)
draw(Hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

# **Univariate Cox Regression Analysis**
library(survival)
library(survminer)

# Define clinical features to be used for univariate Cox regression
Clis <- c("Age", "Gender", "Race", "Child", "Grade", "Stage", "Hepatitis", "Alcohol", "AFP", "PT", "ALB", 'Cluster')
Uni_sur <- sapply(Clis, function(x) as.formula(paste('Surv(Time, Status)~', x)))
Uni_cox <- lapply(Uni_sur, function(x) {coxph(x, data = Cli_data)})

# Extract the results of univariate Cox regression analysis
Uni_results <- lapply(Uni_cox, function(x) {
  x <- summary(x)
  p.value <- signif(x$wald["pvalue"], digits = 2)
  HR <- signif(x$coef[2], digits = 2)
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], digits = 2)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"], digits = 2)
  res <- c(p.value, HR, HR.confint.lower, HR.confint.upper)
  names(res) <- c("p.value", 'HR', 'HR.confint.lower', 'HR.confint.upper')
  return(res)
})

# Combine the results into a data frame
TCGA_Res_uni_cox <- as.data.frame(t(as.data.frame(Uni_results, check.names = FALSE))) %>%
  mutate(ID = paste0('TCGA_', rownames(.)), id = rownames(.))

# **Plot Forest Plot for Univariate Cox Results**
library(ggplot2)
p1 <- ggplot(All_Uni_res, aes(color = P.HR.label, x = HR, xmin = HR.confint.lower, xmax = HR.confint.upper, y = ID)) +
  geom_pointrange(size = 0, fatten = 0) +
  scale_color_manual(values = c("#ce5f06", "#95c9e5", "#0767a9", "#febc84")) +
  theme_bw() +
  coord_cartesian(xlim = c(-1, 6)) +
  xlab("HR (95% CI)") +
  ylab("") +
  theme(legend.position = "none", axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10)) +
  geom_vline(xintercept = 1, color = "black", size = 0.4)

# Save the forest plot as a PDF
ggsave("F:/2023/MOVICS_HCC/Results/All_Uni_forest_plot.pdf", p1, width = 4, height = 8)

# **Plot the Legend for the Forest Plot**
colors <- c("#ce5f06", "#95c9e5", "#0767a9", "#febc84")
names <- c("Positive effect (P<0.05)", "Positive effect (P≥0.05)", "Negative effect (P<0.05)", "Negative effect (P≥0.05)")
pdf("F:/2023/MOVICS_HCC/Results/All_Uni_forest.legend.pdf", width = 6.5, height = 6.0)
plot.new()
legend("top", legend = names, col = colors, pch = 19, horiz = F, pt.cex = 2, box.lwd = 0, box.col = "transparent")
dev.off()
