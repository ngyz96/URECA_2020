library(WGCNA)
library(RColorBrewer)
library(preprocessCore)

bm_10 <- readRDS('Data/02_vst_10_filtered_limma_corrected.rds')
bm_10_coldata <- readRDS('Data/02_bm_10_data.rds')

bm_10_NP.1_coldata <- bm_10_coldata[bm_10_coldata$NP.1 == 1, ]
bm_10_NP.4_coldata <- bm_10_coldata[bm_10_coldata$NP.1 == 4, ]

bm_10_NP.1 <- t(bm_10[, bm_10_NP.1_coldata$index_id]) #rows -> samples, col -> genes
bm_10_NP.4 <- t(bm_10[, bm_10_NP.4_coldata$index_id]) #rows -> samples, col -> genes
    
beta1=3
AdjMatC1<-sign(cor(bm_10_NP.1,method="spearman"))*(cor(bm_10_NP.1,method="spearman"))^2
AdjMatC2<-sign(cor(bm_10_NP.4,method="spearman"))*(cor(bm_10_NP.4,method="spearman"))^2
diag(AdjMatC1)<-0
diag(AdjMatC2)<-0
collectGarbage()

dissTOMC1C2=TOMdist((abs(AdjMatC1-AdjMatC2)/2)^(beta1/2))
collectGarbage()

#Hierarchical clustering is performed using the Topological Overlap of the adjacency difference as input distance matrix
geneTreeC1C2 = flashClust(as.dist(dissTOMC1C2), method = "average");

# Plot the resulting clustering tree (dendrogram)
png(file="hierarchicalTree.png",height=1000,width=1000)
plot(geneTreeC1C2, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);
dev.off()

#We now extract modules from the hierarchical tree. This is done using cutreeDynamic. Please refer to WGCNA package documentation for details
dynamicModsHybridC1C2 = cutreeDynamic(dendro = geneTreeC1C2, distM = dissTOMC1C2,method="hybrid",cutHeight=.996,deepSplit = T, pamRespectsDendro = FALSE,minClusterSize = 20);

#Every module is assigned a color. Note that GREY is reserved for genes which do not belong to any differentially coexpressed module
dynamicColorsHybridC1C2 = labels2colors(dynamicModsHybridC1C2)

#the next step merges clusters which are close (see WGCNA package documentation)
mergedColorC1C2<-mergeCloseModules(rbind(datC1,datC2),dynamicColorsHybridC1C2,cutHeight=.2)$color
colorh1C1C2<-mergedColorC1C2

#reassign better colors
colorh1C1C2[which(colorh1C1C2 =="midnightblue")]<-"red"
colorh1C1C2[which(colorh1C1C2 =="lightgreen")]<-"yellow"
colorh1C1C2[which(colorh1C1C2 =="cyan")]<-"orange"
colorh1C1C2[which(colorh1C1C2 =="lightcyan")]<-"green"
# Plot the dendrogram and colors underneath
png(file="module_assignment.png",width=1000,height=1000)
plotDendroAndColors(geneTreeC1C2, colorh1C1C2, "Hybrid Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors cells")
dev.off()

#We write each module to an individual file containing affymetrix probeset IDs
modulesC1C2Merged<-extractModules(colorh1C1C2,datC1,anno,dir="modules",file_prefix=paste("Output","Specific_module",sep=''),write=T)
write.table(colorh1C1C2,file="module_assignment.txt",row.names=F,col.names=F,quote=F)

#We plot to a file the comparative heatmap showing correlation changes in the modules
#The code for the function plotC1C2Heatmap and others can be found below under the Supporting Functions section
plotC1C2Heatmap(colorh1C1C2,AdjMatC1,AdjMatC2, datC1, datC2)
png(file="exprChange.png",height=500,width=500)
plotExprChange(datC1,datC2,colorh1C1C2)
dev.off()