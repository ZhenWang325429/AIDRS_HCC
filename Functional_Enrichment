########Functional_Enrichment #################

# convert beta value to M value for stronger signal
indata <- LIHC_data[1:4]
indata$meth <- log2(indata$meth / (1 - indata$meth))

# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,F)) # no scale for mutation

feat   <- iClusterBayes.res$feat.res
feat1  <- feat[which(feat$dataset == "mRNA"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "meth"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "mut"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4)

# set color for each omics data
mRNA.col   <- c("#31B29E", "black","#EB6C5A")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
meth.col   <- c("#0099CC", "white","#CC0033")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, meth.col, mut.col)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.TPM","lncRNA.TPM","M value","Mutated"),
             clust.res     = cmoic.LIHC$clust.res, # consensusMOIC Results_V4
             clust.dend    = NULL, # show no dendrogram for samples
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(F,F,F,F), # show no dendrogram for features
             annRow        = annRow, # no selected features
             color         = col.list,
             fig.path      = 'F:/2023/MOVICS_HCC/Results/',
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF CONSENSUSMOIC")

# run GSEA to identify up-regulated GO pathways using Results_V4 from edgeR
gsea.up <- runGSEA(moic.res     = cmoic.LIHC,
                   dea.method   = "edger", # name of DEA method
                   prefix       = "TCGA-LIHC", # MUST be the same of argument in runDEA()
                   dat.path     = 'F:/2023/MOVICS_HCC/Results/', # path of DEA files
                   res.path     = 'F:/2023/MOVICS_HCC/Results/', # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = mRNA_expr_fpkm_symbol_max, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "UPREGULATED PATHWAY HEATMAP")

