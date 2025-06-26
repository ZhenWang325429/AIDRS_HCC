####EaSIeR########

# Load required libraries
library(tidyverse)
library(easier)
library(easierData)

# Custom function to process the data and predict immune responses
EaSIeR_data <- function(Exp_Counts = Exp_Counts, Exp_TPM = Exp_TPM, ID = ID) {
  # Reannotate genes based on provided gene symbols
  genes_info <- easier:::reannotate_genes(cur_genes = rownames(Exp_TPM))
  
  # Remove unsupported gene symbols and entries that are withdrawn
  non_na <- !is.na(genes_info$new_names)
  Exp_TPM <- Exp_TPM[non_na, ]
  genes_info <- genes_info[non_na, ]
  Exp_TPM <- Exp_TPM[!(genes_info$new_names == "entry withdrawn"), ]
  genes_info <- genes_info[!(genes_info$new_names == "entry withdrawn"), ]
  
  # Identify duplicate genes and aggregate them
  newnames_dup <- unique(genes_info$new_names[duplicated(genes_info$new_names)])
  newnames_dup_ind <- do.call(c, lapply(newnames_dup, function(X) which(genes_info$new_names == X)))
  newnames_dup <- genes_info$new_names[newnames_dup_ind]
  
  if (is.null(newnames_dup_ind)) {
    Exp_TPM <- Exp_TPM
  } else {
    tmp <- Exp_TPM[genes_info$old_names[genes_info$new_names %in% newnames_dup], ]
    Exp_TPM <- Exp_TPM[!(rownames(Exp_TPM) %in% rownames(tmp)), ]
    dup_genes <- genes_info$new_names[genes_info$new_names %in% newnames_dup]
    names(dup_genes) <- rownames(tmp)
    tmp2 <- stats::aggregate(tmp, by = list(dup_genes), FUN = "mean")
    rownames(tmp2) <- tmp2$Group.1
    tmp2$Group.1 <- NULL
    Exp_TPM <- rbind(Exp_TPM, tmp2)
  }
  
  # Remove duplicate column names
  Exp_TPM <- Exp_TPM[, !duplicated(colnames(Exp_TPM))]
  Exp_Counts <- Exp_Counts[, !duplicated(colnames(Exp_TPM))]
  
  # Save processed data
  save(Exp_TPM, Exp_Counts, file = paste0('F:/2023/MOVICS_HCC/Data_V4/EaSleR/', ID, '_EaSIeR_data.Rdata'))
  print(paste0(ID, '_Data OK'))
}

# Function to compute immune response scores
MyEaSIeR <- function(Exp_Counts = Exp_Counts, Exp_TPM = Exp_TPM, Cancer = Cancer, ID = ID) {
  # Compute hallmark immune response scores
  hallmarks_of_immune_response <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
  immune_response_scores <- compute_scores_immune_response(RNA_tpm = Exp_TPM, selected_scores = hallmarks_of_immune_response)
  
  # Compute TME features
  cell_fractions <- compute_cell_fractions(RNA_tpm = Exp_TPM)
  pathway_activities <- compute_pathway_activity(RNA_counts = Exp_Counts, remove_sig_genes_immune_response = TRUE)
  
  # Compute TF activities and ligand-receptor pairs
  tf_activities <- compute_TF_activity(RNA_tpm = Exp_TPM)
  lrpair_weights <- compute_LR_pairs(RNA_tpm = Exp_TPM, cancer_type = "pancan")
  ccpair_scores <- compute_CC_pairs(lrpairs = lrpair_weights, cancer_type = "pancan")
  
  # Predict immune response
  predictions <- predict_immune_response(pathways = pathway_activities,
                                         immunecells = cell_fractions,
                                         tfs = tf_activities,
                                         lrpairs = lrpair_weights,
                                         ccpairs = ccpair_scores,
                                         cancer_type = Cancer, 
                                         verbose = TRUE)
  
  # Save the results
  save(hallmarks_of_immune_response, immune_response_scores, pathway_activities, 
       cell_fractions, tf_activities, lrpair_weights, ccpair_scores, predictions, 
       file = paste0('F:/2023/MOVICS_HCC/Data_V4/EaSleR/', ID, '_EaSIeR.Rdata'))
}

# List of IDs to process
IDs = c('GSE144269','GSE141198','GSE141200','GSE124751',
        'ICGC','GSE14520','GSE36376','GSE63898','TCGA_LIHC') 

Counts_list <- list('GSE144269' = Exp_GSE144269_Counts,
                    'GSE141198' = Exp_GSE141198_Counts_filter,
                    'GSE141200' = Exp_GSE141200_Counts_filter,
                    'GSE124751' = Exp_GSE124751_Counts_filter,
                    'ICGC' = Exp_ICGC_Counts_filter,
                    'GSE14520' = Exp_GSE14520_Array_filter,
                    'GSE36376' = Exp_GSE36376_Array_filter,
                    'GSE63898' = Exp_GSE63898_Array_filter,
                    'TCGA_LIHC' = LIHC_Counts)

TPM_list <- list('GSE144269' = Exp_GSE144269_TPM,
                 'GSE141198' = Exp_GSE141198_TPM_filter,
                 'GSE141200' = Exp_GSE141200_TPM_filter,
                 'GSE124751' = Exp_GSE124751_TPM_filter,
                 'ICGC' = Exp_ICGC_TPM_filter,
                 'GSE14520' = Exp_GSE14520_Array_filter,
                 'GSE36376' = Exp_GSE36376_Array_filter,
                 'GSE63898' = Exp_GSE63898_Array_filter,
                 'TCGA_LIHC' = LIHC_TPM)

# Process each dataset
for (i in 1:length(Counts_list)) {
  EaSIeR_data(Exp_Counts = Counts_list[[i]], Exp_TPM = TPM_list[[i]], ID = IDs[i])
}

# Predict immune response for each dataset
for (i in 1:length(IDs)) {
  tryCatch({
    load(paste0('F:/2023/MOVICS_HCC/Data_V4/EaSleR/', IDs[i], '_EaSIeR_data.Rdata'))
    MyEaSIeR(Exp_Counts = Exp_Counts, Exp_TPM = Exp_TPM, Cancer = 'LIHC', ID = IDs[i])
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
  })
}
