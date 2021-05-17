library(WGCNA)

options(stringsAsFactors = FALSE)
vst_counts <- readRDS('Data/02_vst_44_filtered_limma_corrected.rds')
vst_counts <- t(vst_counts)

gsg <- goodSamplesGenes(vst_counts, verbose = 3)
if (!gsg$allOK)
{
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(names(vst_counts)[!gsg$goodGenes], collapse = ", ")))
    if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", paste(rownames(vst_counts)[!gsg$goodSamples], collapse = ", ")))
    # Remove the offending genes and samples from the data:
    vst_counts <- vst_counts[gsg$goodSamples, gsg$goodGenes]
}

sampleTree <- hclust(dist(vst_counts), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "Plots/bm_44_wgcna_clustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

#clinical data
bm_10_coldata <- readRDS('Data/02_bm_44_data.rds')
reduced_coldata <- bm_10_coldata[, c('index_id', 'batch', 'RIN', 'PMI', 'AOD', 'SEX', 'NP.1', 'CDR', 'PlaqueMean', 'bbscore', 'Apo1', 'Apo2')]
reduced_coldata$SEX <- as.factor(reduced_coldata$SEX) #F is 1, M is 2
reduced_coldata$batch <- droplevels(reduced_coldata$batch)
reduced_coldata$batch <- plyr::mapvalues(reduced_coldata$batch, 
                                         from = c('B18C014', 'B28C849', 'B82C014', 'E008C189', 'E009C189', 'K79C014', 'K80C014', 'K82C014', 'K85C014', 'L43C014'), 
                                         to = c(1, 2, 3, 4, 5, 6, 7, 8, 9,10))
reduced_coldata$NP.1 <- plyr::mapvalues(reduced_coldata$NP.1, c('1', '2', '3', '4'), c("1", "4", "3", '2'))
samples <- rownames(vst_counts)
traitRows <- match(samples, reduced_coldata$index_id)
datTraits <- reduced_coldata[traitRows, -1]
rownames(datTraits) <- reduced_coldata[traitRows, 1]
collectGarbage()
save(datTraits, vst_counts, file = "Data/bm_44_wgcna_data_input.RData")

enableWGCNAThreads()
# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft <- pickSoftThreshold(vst_counts, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", 
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(vst_counts, power = 4, maxBlockSize = 10000,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "bm_44TOM",
                       verbose = 3)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "Data/bm_44_wgcna_network_signed.RData")

# Define numbers of genes and samples
nGenes = ncol(vst_counts)
nSamples = nrow(vst_counts)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(vst_counts, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))