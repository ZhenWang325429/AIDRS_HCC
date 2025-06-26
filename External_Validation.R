######## External_Validation #################

# validate with external data using NTP (Yau cohort)
NTP.pred_ICGC <- runNTP(expr = log2(Exp_ICGC_FPKM_filter+1),
                        templates  = marker.up$templates, # the template has been already prepared in runMarker()
                        scaleFlag  = TRUE, # scale input data (by default)
                        centerFlag = TRUE, # center input data (by default)
                        doPlot     = TRUE, # to generate heatmap
                        fig.path   = 'F:/2023/MOVICS_HCC/Results/',
                        fig.name   = "NTP HEATMAP FOR ICGC")

# compare survival outcome in Yau cohort
surv.ICGC <- compSurv(moic.res         = NTP.pred_ICGC ,
                      surv.info        = ICGC_clin.info,
                      convt.time       = "m", # switch to month
                      surv.median.line = "hv", # switch to both
                      fig.name         = "OS KAPLAN-MEIER CURVE OF NTP FOR ICGC")

# You can repeat similar external validation for other datasets like GSE124751, GSE141198, etc.
# The structure is similar: runNTP followed by compSurv.
