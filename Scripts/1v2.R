library(DESeq2)
raw_exp_bm_36 <- read.delim('Data/AMP-AD_MSBB_MSSM_BM_22.raw_counts.tsv') #261
data <- readRDS('Data/01_covariates_cleaned.rds')
samples <- colnames(raw_exp_bm_36)
#keep only caucasians
data <- data[data$RACE == 'W',]
#keep only np.1 = 1 or 4
# data <- data[(data$NP.1 == 1 | data$NP.1 == 2),]

data_bm_36 <- data[data$index_id %in% samples,]
remove_col <- setdiff(colnames(raw_exp_bm_36), data_bm_36$index_id)
raw_exp_bm_36 <- raw_exp_bm_36[, !(colnames(raw_exp_bm_36) %in% remove_col)]
data_bm_36 <- data_bm_36[match(colnames(raw_exp_bm_36), data_bm_36$index_id),]

data_bm_36$AOD <- as.numeric(data_bm_36$AOD)
data_bm_36$AOD <- cut(data_bm_36$AOD, breaks = c(60,69,79,89), labels = c("60-69", "70-79", "80-89"))
levels(data_bm_36$AOD) <- c("60-69", "70-79", "80-89", "90+")
data_bm_36$AOD[is.na(data_bm_36$AOD)] <- "90+"
#binning AOD
data_bm_36$RIN <- cut(data_bm_36$RIN, breaks = 0:10)
#clean data for DESeq object
data_bm_36$NP.1 <- as.factor(data_bm_36$NP.1)
data_bm_36$SEX <- as.factor(data_bm_36$SEX)
data_bm_36$batch <- droplevels(data_bm_36$batch)
#create DESeq object
bm_36 <- DESeqDataSetFromMatrix(countData = raw_exp_bm_36, colData = data_bm_36,
design = ~ NP.1)
keep <- rowSums(counts(bm_36)) >= 200
bm_36 <- bm_36[keep,]
keep <- rowSums(counts(bm_36) > 0) > ncol(counts(bm_36))/2
bm_36 <- bm_36[keep,]
# bm_22_vst <- vst(bm_22)
# plotPCA(bm_22_vst, intgroup = "SEX")
# plotPCA(bm_22_vst, intgroup = "RIN")
# plotPCA(bm_22_vst, intgroup = "NP.1")
# plotPCA(bm_22_vst, intgroup = "batch")
# plotPCA(bm_22_vst, intgroup = "PMI")
# plotPCA(bm_22_vst, intgroup = "AOD")

library("BiocParallel")
register(SnowParam(4))
bm_36 <- DESeq(bm_36)
res <- results(bm_36)
summary(res)
