#load in clinical data
clinical_data <- readxl::read_xlsx('Data/MSBB_clinical.xlsx')
#load in covariates 
covariates <- read.csv('Data/MSBB_RNAseq_covariates_November2018Update.csv') #2559 rows x 16 cols

#filtering out non-okay
covariates <- covariates[covariates$Action == 'OKay',] #2449 rows x 16 cols

#filter out repeated due to bam/fastq files, keeping bam then fastq
covariates <- covariates[order(covariates$fileType),]
covariates <- covariates[!duplicated(covariates[,'sampleIdentifier']),] #1228 rows x 16 cols

#remove '_resequenced' from sampleIdentifier for index creation
index <- grep('resequenced$', covariates$sampleIdentifier)
covariates[index, 'sampleIdentifier'] <- substr(covariates[index, 'sampleIdentifier'],1, nchar(as.character(covariates[index, 'sampleIdentifier'])) - 12) #1228 rows x 16 cols

#remove batch id from end of sampleIdentifier if present
index <- grep('[A-z0-9]{1,}_[A-z0-9]{1,}_[A-z0-9]{1,}_', covariates$sampleIdentifier) #54 samples
covariates[index, 'sampleIdentifier'] <- sub('_[A-Z0-9]{6,}$', '', covariates[index, 'sampleIdentifier'])

#constructing index id for matching to gene matrix
covariates$index_id <- paste0(covariates$batch, '.' ,covariates$sampleIdentifier) #1228 rows x 17 cols

#combine with clinical data
covariates <- merge(covariates, clinical_data, by = 'individualIdentifier') #1228 rows x 27 cols

#check for mismatch data and ensure no missing data after cleaning
length(unlist(lapply(list(raw_exp_bm_10, raw_exp_bm_22, raw_exp_bm_36, raw_exp_bm_44), function(x){setdiff(colnames(x), covariates$index_id)})))
covariates <- covariates[!duplicated(covariates$index_id),] #1213 rows x 27 cols

#remove data individuals with wrong sex
mismatch <- covariates[covariates$SEX.inferred != covariates$SEX,] #14 rows x 27 cols
remove_individuals <- unique(mismatch$individualIdentifier)
covariates <- covariates[!(covariates$individualIdentifier %in% remove_individuals),] #1213 rows x 27 cols

#remove wrongly annotated individual
covariates_m <- covariates[covariates$SEX == 'M',]
raw_exp_bm_10 <- read.delim('Data/AMP-AD_MSBB_MSSM_BM_10.raw_counts.tsv') #261
raw_exp_bm_22 <- read.delim('Data/AMP-AD_MSBB_MSSM_BM_22.raw_counts.tsv') #240
raw_exp_bm_36 <- read.delim('Data/AMP-AD_MSBB_MSSM_BM_36.raw_counts.tsv') #215
raw_exp_bm_44 <- read.delim('Data/AMP-AD_MSBB_MSSM_BM_44.raw_counts.tsv') #222
raw_combined <- cbind(raw_exp_bm_10, raw_exp_bm_22, raw_exp_bm_36, raw_exp_bm_44)
raw_combined_male <- raw_combined[, colnames(raw_combined) %in% covariates_m$index_id]
raw_combined_male_xist <- t(raw_combined_male['ENSG00000229807',]) #XIST gene
hist(raw_combined_male_xist, labels = T)
remove_individuals <- rownames(raw_combined_male_xist)[which.max(raw_combined_male_xist)]
covariates <- covariates[!(covariates$index_id == remove_individual),] #1212 rows x 27 cols

#save the cleaned data
saveRDS(covariates, "Data/covariates_clinical_cleaned.rds")
