###TCIA#########

# Function to calculate Immunophenoscore (IPS)
MyIPS <- function(ID = ID) {
  # Mapping IPS score from expression data
  ipsmap <- function(x) {
    if (x <= 0) {
      ips <- 0
    } else if (x >= 3) {
      ips <- 10
    } else {
      ips <- round(x * 10 / 3, digits = 0)
    }
    return(ips)
  }
  
  # Load expression data for the specified dataset
  load(file = paste0('F:/2023/MOVICS_HCC/Data_V4/EaSleR/', ID, '_EaSIeR_data.Rdata'))
  
  # Apply log2 transformation to TPM data
  if (ID %in% c('GSE144269', 'GSE141198', 'GSE141200', 'GSE124751', 'ICGC', 'TCGA_LIHC')) {
    Exp = apply(Exp_TPM, 2, FUN = function(x) { log2(x + 1) })
  } else {
    Exp = as.matrix(Exp_TPM)
  }
  
  # Read IPS genes and their weights
  IPSG <- read.table(file = 'F:/2023/MOVICS_HCC/Data_V4/EaSleR/Immunophenogram/IPS_genes.txt', header = TRUE, sep = "\t", dec = ".", check.names = FALSE)
  
  # Compute IPS for each sample
  IPS <- NULL
  for (i in 1:length(sample_names)) {
    GE <- Exp[, i]
    mGE <- mean(GE)
    sGE <- sd(GE)
    Z1 <- (as.data.frame(gene_expression)[as.vector(IPSG$GENE), i] - mGE) / sGE
    W1 <- IPSG$WEIGHT
    # Perform calculations for MHC, CP, EC, SC, and final IPS
    # Output results
    # Save final IPS score for each sample
    DF <- data.frame(ID = sample_names, MHC = MHC, EC = EC, SC = SC, CP = CP, AZ = AZ, IPS = IPS)
    write.table(DF, file = paste0('F:/2023/MOVICS_HCC/Data_V4/EaSleR/', ID, "_IPS.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
  }
}

# Run IPS calculations for each dataset
IDs = c('GSE144269', 'GSE141198', 'GSE141200', 'GSE124751', 'ICGC', 'GSE14520', 'GSE36376', 'GSE63898', 'TCGA_LIHC')
for (i in 1:length(IDs)) {
  ID = IDs[i]
  MyIPS(ID)
}
