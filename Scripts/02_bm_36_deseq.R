library(DESeq2)

raw_exp_bm_36 <- read.delim('Data/AMP-AD_MSBB_MSSM_BM_36.raw_counts.tsv') #261
data <- readRDS('Data/01_covariates_cleaned.rds')
samples <- colnames(raw_exp_bm_36)

#keep only caucasians
data <- data[data$RACE == 'W',]

#trimming both covariates data and exp data, keeping only those in common
data_bm_36 <- data[data$index_id %in% samples,]
remove_col <- setdiff(colnames(raw_exp_bm_36), data_bm_36$index_id)
raw_exp_bm_36 <- raw_exp_bm_36[, !(colnames(raw_exp_bm_36) %in% remove_col)]
data_bm_36 <- data_bm_36[match(colnames(raw_exp_bm_36), data_bm_36$index_id),]
saveRDS(data_bm_36, "Data/02_bm_36_data.rds")

#create matrix for WGCNA
bm_36_wgcna <- DESeqDataSetFromMatrix(countData = raw_exp_bm_36, colData = data_bm_36, 
                                      design = ~batch + RIN + SEX + PMI + AOD + NP.1)
keep <- rowSums(counts(bm_36_wgcna)) >= 200
bm_36_wgcna <- bm_36_wgcna[keep,]
keep <- rowSums(counts(bm_36_wgcna) > 0) > ncol(counts(bm_36_wgcna))/2
bm_36_wgcna <- bm_36_wgcna[keep,]
mat <- assay(vst(bm_36_wgcna))
mat <- limma::removeBatchEffect(mat, batch = bm_36_wgcna$batch, batch2 = bm_36_wgcna$SEX, covariates = cbind(bm_36_wgcna$RIN,bm_36_wgcna$PMI,bm_36_wgcna$AOD))
saveRDS(mat, file = 'Data/02_vst_36_filtered_limma_corrected.rds')

#binning age groups
data_bm_36$AOD <- as.numeric(data_bm_36$AOD)
data_bm_36$AOD <- cut(data_bm_36$AOD, breaks = c(60,69,79,89), labels = c("60-69", "70-79", "80-89"))
levels(data_bm_36$AOD) <- c("60-69", "70-79", "80-89", "90+")
data_bm_36$AOD[is.na(data_bm_36$AOD)] <- "90+"

#binning AOD
data_bm_36$RIN <- cut(data_bm_36$RIN, breaks = 0:10)

# #binning PMI
# data$PMI <- cut(data$PMI, breaks = c(0,200,400,600,800,1000,1200,1400,1600,1800),
#                 labels = c("0-200","200-400", "400-600","600-800","800-1000","1000-1200","1200-1400","1400-1600","1600-1800"))

#clean data for DESeq object
data_bm_36$NP.1 <- as.factor(data_bm_36$NP.1)
data_bm_36$SEX <- as.factor(data_bm_36$SEX)
data_bm_36$batch <- droplevels(data_bm_36$batch)

#create DESeq object
bm_36 <- DESeqDataSetFromMatrix(countData = raw_exp_bm_36, colData = data_bm_36, 
                                design = ~batch + RIN + SEX + PMI + AOD + NP.1)
keep <- rowSums(counts(bm_36)) >= 200
bm_36 <- bm_36[keep,]
keep <- rowSums(counts(bm_36) > 0) > ncol(counts(bm_36))/2
bm_36 <- bm_36[keep,]
bm_36_vst <- vst(bm_36)
plotPCA(bm_36_vst, intgroup = "SEX")
saveRDS(assay(bm_36_vst), file = 'Data/02_vst_36_filtered.rds')

#visualise all design factors
plotPCA(bm_36_vst, intgroup = "SEX")
plotPCA(bm_36_vst, intgroup = "RIN")
plotPCA(bm_36_vst, intgroup = "NP.1")
plotPCA(bm_36_vst, intgroup = "batch")
plotPCA(bm_36_vst, intgroup = "PMI")
plotPCA(bm_36_vst, intgroup = "AOD")

#bm_36 <- DESeq(bm_36, test = 'LRT', reduced = ~batch + RIN + SEX + PMI + AOD)
bm_36 <- DESeq(bm_36)
bm_36_result <- results(bm_36, name="NP.1_2_vs_1")
summary(bm_36_result)
saveRDS(bm_36, 'Data/02_bm_36_Wald.rds')
