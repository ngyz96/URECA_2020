library(DESeq2)

raw_exp_bm_22 <- read.delim('Data/AMP-AD_MSBB_MSSM_BM_22.raw_counts.tsv') #261
data <- readRDS('Data/01_covariates_cleaned.rds')
samples <- colnames(raw_exp_bm_22)

#keep only caucasians
data <- data[data$RACE == 'W',]

#trimming both covariates data and exp data, keeping only those in common
data_bm_22 <- data[data$index_id %in% samples,]
remove_col <- setdiff(colnames(raw_exp_bm_22), data_bm_22$index_id)
raw_exp_bm_22 <- raw_exp_bm_22[, !(colnames(raw_exp_bm_22) %in% remove_col)]
data_bm_22 <- data_bm_22[match(colnames(raw_exp_bm_22), data_bm_22$index_id),]
saveRDS(data_bm_22, "Data/02_bm_22_data.rds")

#create matrix for WGCNA
bm_22_wgcna <- DESeqDataSetFromMatrix(countData = raw_exp_bm_22, colData = data_bm_22, 
                                      design = ~batch + RIN + SEX + PMI + AOD + NP.1)
keep <- rowSums(counts(bm_22_wgcna)) >= 200
bm_22_wgcna <- bm_22_wgcna[keep,]
keep <- rowSums(counts(bm_22_wgcna) > 0) > ncol(counts(bm_22_wgcna))/2
bm_22_wgcna <- bm_22_wgcna[keep,]
mat <- assay(vst(bm_22_wgcna))
mat <- limma::removeBatchEffect(mat, batch = bm_22_wgcna$batch, batch2 = bm_22_wgcna$SEX, covariates = cbind(bm_22_wgcna$RIN,bm_22_wgcna$PMI,bm_22_wgcna$AOD))
saveRDS(mat, file = 'Data/02_vst_22_filtered_limma_corrected.rds')

#binning age groups
data_bm_22$AOD <- as.numeric(data_bm_22$AOD)
data_bm_22$AOD <- cut(data_bm_22$AOD, breaks = c(60,69,79,89), labels = c("60-69", "70-79", "80-89"))
levels(data_bm_22$AOD) <- c("60-69", "70-79", "80-89", "90+")
data_bm_22$AOD[is.na(data_bm_22$AOD)] <- "90+"

# #binning PMI
# data$PMI <- cut(data$PMI, breaks = c(0,200,400,600,800,1000,1200,1400,1600,1800),
#                 labels = c("0-200","200-400", "400-600","600-800","800-1000","1000-1200","1200-1400","1400-1600","1600-1800"))

#binning AOD
data_bm_22$RIN <- cut(data_bm_22$RIN, breaks = 0:10)

#clean data for DESeq object
data_bm_22$NP.1 <- as.factor(data_bm_22$NP.1)
data_bm_22$SEX <- as.factor(data_bm_22$SEX)
data_bm_22$batch <- droplevels(data_bm_22$batch)

#create DESeq object
bm_22 <- DESeqDataSetFromMatrix(countData = raw_exp_bm_22, colData = data_bm_22, 
                                design = ~batch + RIN + SEX + PMI + AOD + NP.1)
keep <- rowSums(counts(bm_22)) >= 200
bm_22 <- bm_22[keep,]
keep <- rowSums(counts(bm_22) > 0) > ncol(counts(bm_22))/2
bm_22 <- bm_22[keep,]
bm_22_vst <- vst(bm_22)
plotPCA(bm_22_vst, intgroup = "SEX")
saveRDS(assay(bm_22_vst), file = 'Data/02_vst_22_filtered.rds')

#visualise all design factors
plotPCA(bm_22_vst, intgroup = "SEX")
plotPCA(bm_22_vst, intgroup = "RIN")
plotPCA(bm_22_vst, intgroup = "NP.1")
plotPCA(bm_22_vst, intgroup = "batch")
plotPCA(bm_22_vst, intgroup = "PMI")
plotPCA(bm_22_vst, intgroup = "AOD")

#bm_22 <- DESeq(bm_22, test = 'LRT', reduced = ~batch + RIN + SEX + PMI + AOD)
bm_22 <- DESeq(bm_22)
bm_22_result <- results(bm_22, name="NP.1_2_vs_1")
summary(bm_22_result)
saveRDS(bm_22, 'Data/02_bm_22_WALD.rds')
