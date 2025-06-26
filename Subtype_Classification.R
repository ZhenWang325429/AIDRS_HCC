######## 1. Subtype #################

# identify optimal clustering number (may take a while)
optk.LIHC <- getClustNum(data        = LIHC_data[1:4],
                         is.binary   = c(F,F,F,T), # note: the 4th data is somatic mutation which is a binary matrix
                         try.N.clust = 2:8, # try cluster number from 2 to 8
                         fig.name    = "CLUSTER NUMBER OF TCGA-LIHC",
                         fig.path    = 'F:/2023/MOVICS_HCC/Results/')

# **分型**
iClusterBayes.res <- getMOIC(data        = LIHC_data[1:4],
                             N.clust     = 2,
                             methodslist = "iClusterBayes", # specify only ONE algorithm here
                             type        = c("gaussian","gaussian","gaussian","binomial"), # data type corresponding to the list
                             n.burnin    = 1800,
                             n.draw      = 1200,
                             prior.gamma = c(0.5, 0.5,  0.5, 0.5),
                             sdev        = 0.05,
                             thin        = 2)

# perform multi-omics integrative clustering with the rest of 9 algorithms
moic.res.list <- getMOIC(data        = LIHC_data[1:4],
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"),
                         N.clust     = 2,
                         type        = c("gaussian", "gaussian", "gaussian", "binomial"))

# 将第一种算法的分型结果纳入总列表
moic.res.list <- append(moic.res.list, 
                        list("iClusterBayes" = iClusterBayes.res))

# **非常值得注意的是此时的结果相当于是将额外的一种添加进入原有的9种计算结果之中**

cmoic.LIHC <- getConsensusMOIC(moic.res.list = moic.res.list[1:9],
                               fig.path      = 'F:/2023/MOVICS_HCC/Results/',
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")

save(iClusterBayes.res,moic.res.list,cmoic.LIHC,
     file = here('F:/2023/MOVICS_HCC/Results/','04_LIHC_moic.res.list.Rdata'))

getSilhouette(sil      = cmoic.LIHC$sil, # a sil object returned by getConsensusMOIC()
              fig.path = 'F:/2023/MOVICS_HCC/Results/',
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)

