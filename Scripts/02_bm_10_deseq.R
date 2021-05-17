library(DESeq2)
library(BiocParallel)
register(SnowParam(4))

raw_exp_bm_10 <- read.delim('Data/AMP-AD_MSBB_MSSM_BM_10.raw_counts.tsv') #261
data <- readRDS('Data/01_covariates_cleaned.rds')
samples <- colnames(raw_exp_bm_10)

#keep only caucasians
data <- data[data$RACE == 'W',]
#keep only np.1 = 1 or 4
#data <- data[(data$NP.1 == 1 | data$NP.1 == 4),]
#trimming both covariates data and exp data, keeping only those in common
data_bm_10 <- data[data$index_id %in% samples,]
remove_col <- setdiff(colnames(raw_exp_bm_10), data_bm_10$index_id)
raw_exp_bm_10 <- raw_exp_bm_10[, !(colnames(raw_exp_bm_10) %in% remove_col)]
data_bm_10 <- data_bm_10[match(colnames(raw_exp_bm_10), data_bm_10$index_id),]
saveRDS(data_bm_10, "Data/02_bm_10_data.rds")

#create matrix for WGCNA
bm_10_wgcna <- DESeqDataSetFromMatrix(countData = raw_exp_bm_10, colData = data_bm_10, 
                                      design = ~batch + RIN + SEX + PMI + AOD + NP.1)
keep <- rowSums(counts(bm_10_wgcna)) >= 200
bm_10_wgcna <- bm_10_wgcna[keep,]
keep <- rowSums(counts(bm_10_wgcna) > 0) > ncol(counts(bm_10_wgcna))/2
bm_10_wgcna <- bm_10_wgcna[keep,]
bm_10_wgcna_vst <- vst(bm_10_wgcna)
wgcna_pca <- plotPCA(bm_10_wgcna_vst, intgroup = "SEX", returnData = T)
mislabeled <- as.character(wgcna_pca[wgcna_pca$PC1 > 20 & wgcna_pca$SEX == 'M',]$name)
bm_10_wgcna <- bm_10_wgcna[, !(colnames(bm_10_wgcna) == mislabeled)]
mat <- assay(vst(bm_10_wgcna))
mat <- limma::removeBatchEffect(mat, batch = bm_10_wgcna$batch, batch2 = bm_10_wgcna$SEX, covariates = cbind(bm_10_wgcna$RIN,bm_10_wgcna$PMI,bm_10_wgcna$AOD))
saveRDS(mat, file = 'Data/02_vst_10_filtered_limma_corrected.rds')

#binning age groups
data_bm_10$AOD <- as.numeric(data_bm_10$AOD)
data_bm_10$AOD <- cut(data_bm_10$AOD, breaks = c(60,69,79,89), labels = c("60-69", "70-79", "80-89"))
levels(data_bm_10$AOD) <- c("60-69", "70-79", "80-89", "90+")
data_bm_10$AOD[is.na(data_bm_10$AOD)] <- "90+"

# #binning PMI
# data$PMI <- cut(data$PMI, breaks = c(0,200,400,600,800,1000,1200,1400,1600,1800),
#                 labels = c("0-200","200-400", "400-600","600-800","800-1000","1000-1200","1200-1400","1400-1600","1600-1800"))

#binning AOD
data_bm_10$RIN <- cut(data_bm_10$RIN, breaks = 0:10)

#clean data for DESeq object
data_bm_10$NP.1 <- as.factor(data_bm_10$NP.1)
data_bm_10$SEX <- as.factor(data_bm_10$SEX)
data_bm_10$batch <- droplevels(data_bm_10$batch)

#create DESeq object
bm_10 <- DESeqDataSetFromMatrix(countData = raw_exp_bm_10, colData = data_bm_10, 
                                design = ~batch + RIN + SEX + PMI + AOD + NP.1)
keep <- rowSums(counts(bm_10)) >= 200
bm_10 <- bm_10[keep,]
keep <- rowSums(counts(bm_10) > 0) > ncol(counts(bm_10))/2
bm_10 <- bm_10[keep,]
bm_10_vst <- vst(bm_10)
plotPCA(bm_10_vst, intgroup = "SEX")

#remove mislabeled male
pca <- plotPCA(bm_10_vst, intgroup = "SEX", returnData = T)
mislabeled <- as.character(pca[pca$PC1 > 20 & pca$SEX == 'M',]$name)
bm_10 <- bm_10[, !(colnames(bm_10) == mislabeled)]

bm_10_vst <- vst(bm_10)
mat <- assay(bm_10_vst)
saveRDS(assay(bm_10_vst), file = 'Data/02_vst_10_filtered.rds')

#visualise all design factors
plotPCA(bm_10_vst, intgroup = "SEX")
plotPCA(bm_10_vst, intgroup = "RIN")
plotPCA(bm_10_vst, intgroup = "NP.1")
plotPCA(bm_10_vst, intgroup = "batch")
plotPCA(bm_10_vst, intgroup = "PMI")
plotPCA(bm_10_vst, intgroup = "AOD")

# #merge pca data with covariates
# pca_data <- plotPCA(bm_10_vst, intgroup = "batch", returnData = T)
# data_bm_10 <- data_bm_10[!(data_bm_10$index_id == mislabeled), ]
# pca_data <- pca_data[, c('PC1', 'PC2', 'name')]
# colnames(pca_data) <- c('PC1', 'PC2', 'index_id')
# data_with_pca <- merge(data_bm_10, pca_data, by='index_id')

# bm_22 <- DESeq(bm_22, test = 'LRT', reduced = ~batch + RIN + SEX + PMI + AOD)
bm_10 <- DESeq(bm_10)
bm_10_result <- results(bm_10, name="NP.1_2_vs_1")
summary(bm_10_result)
saveRDS(bm_10, 'Data/02_bm_10_LRT.rds')