library(DESeq2)

raw_exp_bm_44 <- read.delim('Data/AMP-AD_MSBB_MSSM_BM_44.raw_counts.tsv') #261
data <- readRDS('Data/01_covariates_cleaned.rds')
samples <- colnames(raw_exp_bm_44)

#keep only caucasians
data <- data[data$RACE == 'W',]

#trimming both covariates data and exp data, keeping only those in common
data_bm_44 <- data[data$index_id %in% samples,]
remove_col <- setdiff(colnames(raw_exp_bm_44), data_bm_44$index_id)
raw_exp_bm_44 <- raw_exp_bm_44[, !(colnames(raw_exp_bm_44) %in% remove_col)]
data_bm_44 <- data_bm_44[match(colnames(raw_exp_bm_44), data_bm_44$index_id),]
saveRDS(data_bm_44, "Data/02_bm_44_data.rds")

#create matrix for WGCNA
bm_44_wgcna <- DESeqDataSetFromMatrix(countData = raw_exp_bm_44, colData = data_bm_44, 
                                design = ~batch + RIN + SEX + PMI + AOD + NP.1)
keep <- rowSums(counts(bm_44_wgcna)) >= 200
bm_44_wgcna <- bm_44_wgcna[keep,]
keep <- rowSums(counts(bm_44_wgcna) > 0) > ncol(counts(bm_44_wgcna))/2
bm_44_wgcna <- bm_44_wgcna[keep,]
mat <- assay(vst(bm_44_wgcna))
mat <- limma::removeBatchEffect(mat, batch = bm_44_wgcna$batch, batch2 = bm_44_wgcna$SEX, covariates = cbind(bm_44_wgcna$RIN,bm_44_wgcna$PMI,bm_44_wgcna$AOD))
saveRDS(mat, file = 'Data/02_vst_44_filtered_limma_corrected.rds')

#binning age groups
data_bm_44$AOD <- as.numeric(data_bm_44$AOD)
data_bm_44$AOD <- cut(data_bm_44$AOD, breaks = c(60,69,79,89), labels = c("60-69", "70-79", "80-89"))
levels(data_bm_44$AOD) <- c("60-69", "70-79", "80-89", "90+")
data_bm_44$AOD[is.na(data_bm_44$AOD)] <- "90+"

#binning AOD
data_bm_44$RIN <- cut(data_bm_44$RIN, breaks = 0:10)

# #binning PMI
# data$PMI <- cut(data$PMI, breaks = c(0,200,400,600,800,1000,1200,1400,1600,1800),
#                 labels = c("0-200","200-400", "400-600","600-800","800-1000","1000-1200","1200-1400","1400-1600","1600-1800"))

#cleaning data for DESeq object
data_bm_44$NP.1 <- as.factor(data_bm_44$NP.1)
data_bm_44$SEX <- as.factor(data_bm_44$SEX)
data_bm_44$batch <- droplevels(data_bm_44$batch)

#create DESeq object
bm_44 <- DESeqDataSetFromMatrix(countData = raw_exp_bm_44, colData = data_bm_44, 
                                design = ~batch + RIN + SEX + PMI + AOD + NP.1)
keep <- rowSums(counts(bm_44)) >= 200
bm_44 <- bm_44[keep,]
keep <- rowSums(counts(bm_44) > 0) > ncol(counts(bm_44))/2
bm_44 <- bm_44[keep,]
bm_44_vst <- vst(bm_44)
plotPCA(bm_44_vst, intgroup = "SEX")
saveRDS(assay(bm_44_vst), file = 'Data/02_vst_44_filtered.rds')

#visualise all design factors
plotPCA(bm_44_vst, intgroup = "SEX")
plotPCA(bm_44_vst, intgroup = "RIN")
plotPCA(bm_44_vst, intgroup = "NP.1")
plotPCA(bm_44_vst, intgroup = "batch")
plotPCA(bm_44_vst, intgroup = "PMI")
plotPCA(bm_44_vst, intgroup = "AOD")

#bm_44 <- DESeq(bm_44, test = 'LRT', reduced = ~batch + RIN + SEX + PMI + AOD)
bm_44 <- DESeq(bm_44)
bm_44_result <- results(bm_44, name="NP.1_2_vs_1")
summary(bm_44_result)
saveRDS(bm_44, 'Data/02_bm_44_Wald.rds')
