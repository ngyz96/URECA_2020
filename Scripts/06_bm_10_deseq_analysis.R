library(gprofiler2)
library(DESeq2)

bm10 <- readRDS("Data/02_bm_10_Wald.rds")
bm22 <- readRDS("Data/02_bm_22_Wald.rds")
bm36 <- readRDS("Data/02_bm_36_Wald.rds")
bm44 <- readRDS("Data/02_bm_44_Wald.rds")

bm10res <- results(bm10, name="NP.1_2_vs_1")
bm22res <- results(bm22, name="NP.1_2_vs_1")
bm36res <- results(bm36, name="NP.1_2_vs_1")
bm44res <- results(bm44, name="NP.1_2_vs_1")
allGenes10 <- rownames(subset(bm10res))
allGenes22 <- rownames(subset(bm22res))
allGenes36 <- rownames(subset(bm36res))
allGenes44 <- rownames(subset(bm44res))
sigGenes10 <- rownames(subset(bm10res, padj < 0.1))
sigGenes22 <- rownames(subset(bm22res, padj < 0.1))
sigGenes36 <- rownames(subset(bm36res, padj < 0.1))
sigGenes44 <- rownames(subset(bm44res, padj < 0.1))

sigGenes36_upregulated <- subset(bm36res, padj < 0.1)
sigGenes36_upregulated <- subset(sigGenes36_upregulated, log2FoldChange > 0)
sigGenes36_upregulated <- sigGenes36_upregulated[order(sigGenes36_upregulated$padj),]
sigGenes36_upregulated_names <- rownames(sigGenes36_upregulated)

sigGenes36_downregulated <- subset(bm36res, padj < 0.1)
sigGenes36_downregulated <- subset(sigGenes36_downregulated, log2FoldChange < 0)
sigGenes36_downregulated <- sigGenes36_downregulated[order(sigGenes36_downregulated$padj),]
sigGenes36_downregulated_names <- rownames(sigGenes36_downregulated)

sigGenes22_upregulated <- subset(bm22res, padj < 0.1)
sigGenes22_upregulated <- subset(sigGenes22_upregulated, log2FoldChange > 0)
sigGenes22_upregulated <- sigGenes22_upregulated[order(sigGenes22_upregulated$padj),]
sigGenes22_upregulated_names <- rownames(sigGenes22_upregulated)

sigGenes22_downregulated <- subset(bm22res, padj < 0.1)
sigGenes22_downregulated <- subset(sigGenes22_downregulated, log2FoldChange < 0)
sigGenes22_downregulated <- sigGenes22_downregulated[order(sigGenes22_downregulated$padj),]
sigGenes22_downregulated_names <- rownames(sigGenes22_downregulated)
gostres22up_regulated <- gost(query = sigGenes22_upregulated_names,
                              organism = "hsapiens", ordered_query = TRUE, 
                              multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                              measure_underrepresentation = FALSE, evcodes = TRUE, 
                              user_threshold = 0.05, correction_method = "g_SCS", 
                              domain_scope = "annotated", custom_bg = NULL, 
                              numeric_ns = "", sources = c("GO", "KEGG", "REAC"))

gostres36up_regulated <- gost(query = sigGenes36_upregulated_names,
                    organism = "hsapiens", ordered_query = TRUE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = TRUE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = c("GO", "KEGG", "REAC"))
gostres36down_regulated <- gost(query = sigGenes36_downregulated_names,
                              organism = "hsapiens", ordered_query = TRUE, 
                              multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                              measure_underrepresentation = FALSE, evcodes = TRUE, 
                              user_threshold = 0.05, correction_method = "g_SCS", 
                              domain_scope = "annotated", custom_bg = NULL, 
                              numeric_ns = "", sources = c("GO", "KEGG", "REAC"))

sigGenes44_upregulated <- subset(bm44res, padj < 0.1)
sigGenes44_upregulated <- subset(sigGenes44_upregulated, log2FoldChange > 0)
sigGenes44_upregulated <- sigGenes44_upregulated[order(sigGenes44_upregulated$padj),]
sigGenes44_upregulated_names <- rownames(sigGenes44_upregulated)

sigGenes44_downregulated <- subset(bm44res, padj < 0.1)
sigGenes44_downregulated <- subset(sigGenes44_downregulated, log2FoldChange < 0)
sigGenes44_downregulated <- sigGenes44_downregulated[order(sigGenes44_downregulated$padj),]
sigGenes44_downregulated_names <- rownames(sigGenes44_downregulated)

gostres44up_regulated <- gost(query = sigGenes44_upregulated_names,
                              organism = "hsapiens", ordered_query = TRUE, 
                              multi_query = FALSE, significant = T, exclude_iea = FALSE, 
                              measure_underrepresentation = FALSE, evcodes = TRUE, 
                              user_threshold = 0.05, correction_method = "g_SCS", 
                              domain_scope = "annotated", custom_bg = NULL, 
                              numeric_ns = "", sources = NULL)
gostres44down_regulated <- gost(query = sigGenes44_downregulated_names,
                                organism = "hsapiens", ordered_query = TRUE, 
                                multi_query = FALSE, significant = T, exclude_iea = FALSE, 
                                measure_underrepresentation = FALSE, evcodes = TRUE, 
                                user_threshold = 0.05, correction_method = "g_SCS", 
                                domain_scope = "annotated", custom_bg = NULL, 
                                numeric_ns = "", sources = c("GO", "KEGG", "REAC"))

gostres10 <- gost(query = sigGenes10,
                 organism = "hsapiens", ordered_query = FALSE, 
                 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                 measure_underrepresentation = FALSE, evcodes = TRUE, 
                 user_threshold = 0.05, correction_method = "g_SCS", 
                 domain_scope = "annotated", custom_bg = NULL, 
                 numeric_ns = "", sources = NULL)
gostres22 <- gost(query = sigGenes22,
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL)
gostres36up <- gost(query = sigGenes36_upregulated_names,
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL)
gostres44 <- gost(query = sigGenes44,
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL)

gostplot(gostres10, capped = TRUE, interactive = TRUE)
gostplot(gostres22, capped = TRUE, interactive = TRUE)
gostplot(gostres36, capped = TRUE, interactive = TRUE)
gostplot(gostres44, capped = TRUE, interactive = TRUE)

gostres36terms <- c('GO:0007268', 'GO:0099536', 'GO:0098916', 'GO:0099537', 'GO:0099003', 
                    'REAC:R-HSA-264642', 'REAC:R-HSA-112316', 'REAC:R-HSA-1428517')
gostres36upterms <- c('GO:0030198', 'GO:0043062', 'GO:0022610', 'GO:0007155', 'GO:0031589',
                      'REAC:R-HSA-1474244', 'REAC:R-HSA-3000178', 'REAC:R-HSA-216083')
gostres44downterms <- c('GO:0022613', 'GO:0032543', 'GO"0034660', 'GO:0042254', 'GO:0006396',
                        'REAC:R-HSA-72766', 'REAC:R-HSA-611105', 'KEGG:00970')
c(gostbmdown$result$term_name[1:5], 'REAC:R-HSA-1428517', 'REAC:R-HSA-163200', 'REAC:R-HSA-611105')