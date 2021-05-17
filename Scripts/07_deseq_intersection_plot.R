library(DESeq2)
library(UpSetR)
library(ggplot2)

bm10 <- readRDS("Data/02_bm_10_Wald.rds")
bm22 <- readRDS("Data/02_bm_22_Wald.rds")
bm36 <- readRDS("Data/02_bm_36_Wald.rds")
bm44 <- readRDS("Data/02_bm_44_Wald.rds")

bm10res <- results(bm10, name="NP.1_2_vs_1")
bm22res <- results(bm22, name="NP.1_2_vs_1")
bm36res <- results(bm36, name="NP.1_2_vs_1")
bm44res <- results(bm44, name="NP.1_2_vs_1")

sigGenes10_upregulated <- rownames(subset(subset(bm10res, padj < 0.1), log2FoldChange > 0))
sigGenes22_upregulated <- rownames(subset(subset(bm22res, padj < 0.1), log2FoldChange > 0))
sigGenes36_upregulated <- rownames(subset(subset(bm36res, padj < 0.1), log2FoldChange > 0))
sigGenes44_upregulated <- rownames(subset(subset(bm44res, padj < 0.1), log2FoldChange > 0))
sigGenes44_downregulated <- rownames(subset(subset(bm44res, padj < 0.1), log2FoldChange < 0))
sigGenes36_downregulated <- rownames(subset(subset(bm36res, padj < 0.1), log2FoldChange < 0))
sigGenes22_downregulated <- rownames(subset(subset(bm22res, padj < 0.1), log2FoldChange < 0))
sigGenes10_downregulated <- rownames(subset(subset(bm10res, padj < 0.1), log2FoldChange < 0))

upregulatedList <- list(BM10 = sigGenes10_upregulated,
                        BM22 = sigGenes22_upregulated,
                        BM36 = sigGenes36_upregulated,
                        BM44 = sigGenes44_upregulated)
downregulatedList <- list(BM10 = sigGenes10_downregulated,
                          BM22 = sigGenes22_downregulated,
                          BM36 = sigGenes36_downregulated,
                          BM44 = sigGenes44_downregulated)
#plotting of upset plots
upset(fromList(upregulatedList), 
      mb.ratio = c(0.65,0.35), 
      mainbar.y.label = "Upregulated Gene Intersection Size",
      sets.x.label = "Number of DEGs",
      point.size = 4.5,
      line.size = 2,
      text.scale = c(3, 1.8, 1.8, 1.8, 2.8, 2.8))
upset(fromList(downregulatedList), 
      mb.ratio = c(0.65,0.35), 
      mainbar.y.label = "Downregulated Gene Intersection Size",
      sets.x.label = "Number of DEGs",
      point.size = 4.5,
      line.size = 2,
      text.scale = c(3, 1.8, 1.8, 1.8, 2.8, 2.8))
