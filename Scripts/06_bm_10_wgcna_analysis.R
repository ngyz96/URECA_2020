library(WGCNA)
options(stringsAsFactors = FALSE)

load(file = "Data/bm_44_wgcna_data_input.RData")
load(file = "Data/bm_44_wgcna_network.RData")

# Define numbers of genes and samples
nGenes = ncol(vst_counts)
nSamples = nrow(vst_counts)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(vst_counts, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Define variable weight containing the weight column of datTrait
CDR <- as.data.frame(datTraits$CDR)
names(CDR) = "CDR"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(vst_counts, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(vst_counts, CDR, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(CDR), sep="")
names(GSPvalue) = paste("p.GS.", names(CDR), sep="")

module = "skyblue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for CDR",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
library(gprofiler2)

module = "darkgreen"
column = match(module, modNames);
moduleGenes = moduleColors==module;
blue<- colnames(vst_counts[,moduleGenes])

gost_blue <- gost(query = blue,
                  organism = "hsapiens", ordered_query = FALSE, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                  measure_underrepresentation = FALSE, evcodes = TRUE, 
                  user_threshold = 0.05, correction_method = "g_SCS", 
                  domain_scope = "annotated", custom_bg = NULL, 
                  numeric_ns = "", sources = c("GO", "KEGG", "REAC"))
gostres44blueterms <- c('KEGG:04020','KEGG:04360', 'REAC:R-HSA-112316', 'REAC:R-HSA-112315', 'REAC:R-HSA-112314') 
