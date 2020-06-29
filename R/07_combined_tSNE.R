#######################################################
# Purpose: Combined t-SNE of RNASeq and Microarray data
# Author: Komal S. Rathi
# Date: 06/28/2020
#######################################################

# load libraries
library(ggplot2)
library(Rtsne)
library(tidyverse)
library(ggpubr)
source('R/utils/pubTheme.R')

# merged GERs
allGeneCombos <- readRDS('data/model/bestFeaturesNew.RDS')

# Figure 2B: T-SNE for DS1
ds1GER <- readRDS("data/RNASeqDataForPlotDS1.RDS")
ds1geneRatioOut <- ds1GER[[1]]
ds1sampAnnot <- ds1GER[[2]]
ds1sampAnnot$Type <- 'RNA-Seq'

# Figure 2C: T-SNE for DS2
ds2GER <- readRDS("data/RNASeqDataForPlotDS2.RDS")
ds2geneRatioOut <- ds2GER[[1]]
ds2sampAnnot <- ds2GER[[2]]
ds2sampAnnot$Type <- 'Microarray'
colnames(ds2sampAnnot) <- colnames(ds1sampAnnot)
ds2sampAnnot$Sample <- rownames(ds2sampAnnot)

# combine dataset
sampAnnot <- rbind(ds1sampAnnot, ds2sampAnnot)
sampAnnot$Subgroup <- as.character(sampAnnot$Subgroup)
sampAnnot$Subgroup[sampAnnot$Subgroup == "U"] <- 'Unknown'
common.genes <- intersect(rownames(ds1geneRatioOut), rownames(ds2geneRatioOut))
common.genes <- intersect(allGeneCombos, common.genes)
ds1geneRatioOut <- ds1geneRatioOut[common.genes,]
ds2geneRatioOut <- ds2geneRatioOut[common.genes,]
mat <- ds1geneRatioOut %>%
  rownames_to_column('ger') %>%
  inner_join(ds2geneRatioOut %>%
               rownames_to_column('ger'), by = 'ger')  %>%
  column_to_rownames('ger')

# 2D t-SNE
set.seed(42)
tsneOut <- Rtsne(t(log2(mat)), initial_dims=200, perplexity=10, max_iter=500)
tsneOut <- data.frame(tsneOut$Y, sampAnnot)
tsneOut$Subgroup <- factor(tsneOut$Subgroup, levels = c("Group3", "Group4", "SHH", "WNT", "Unknown"))
p <- ggplot(tsneOut, aes(X1, X2, shape = Type, color = Subgroup))+
  geom_point(size = 5, alpha = 0.6)+
  theme_bw()+
  theme_Publication(base_size = 10) + xlab("tSNE dimension 1") + ylab("tSNE dimension 2")  +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  scale_color_manual(values = c("Group3" = "#F8766D",
                                "Group4" = "#7CAE00",
                                "SHH" = "#00BFC4",
                                "WNT" = "#C77CFF",
                                "Unknown" = "#000000"))
ggsave(plot = p, device = pdf, filename = "results/plots/2DtSNE_Combined.pdf", width = 7, height = 6)

tsneOut$Group <- paste0(tsneOut$Type, '_',tsneOut$Subgroup)
tsneOut$Group <- factor(tsneOut$Group, levels = c("RNA-Seq_Group3", "RNA-Seq_Group4", "RNA-Seq_SHH", "RNA-Seq_WNT", 
                                                  "Microarray_Group3", "Microarray_Group4", "Microarray_SHH", "Microarray_WNT", "Microarray_Unknown"))
q <- ggplot(tsneOut, aes(X1, X2, shape = Group, color = Group))+
  geom_point(size = 5, alpha = 0.6)+
  theme_bw()+
  theme_Publication(base_size = 10) + xlab("tSNE dimension 1") + ylab("tSNE dimension 2")  +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  scale_colour_manual(name = "Group",
                      labels = c("RNA-Seq_Group3", "RNA-Seq_Group4", "RNA-Seq_SHH", "RNA-Seq_WNT", 
                                 "Microarray_Group3", "Microarray_Group4", "Microarray_SHH", "Microarray_WNT", "Microarray_Unknown"),
                      values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF",
                                 "#F8766D", "#7CAE00", "#00BFC4", "#C77CFF","#000000")) +
  scale_shape_manual(name = "Group",
                     labels = c("RNA-Seq_Group3", "RNA-Seq_Group4", "RNA-Seq_SHH", "RNA-Seq_WNT", 
                                "Microarray_Group3", "Microarray_Group4", "Microarray_SHH", "Microarray_WNT", "Microarray_Unknown"),
                     values = c(17, 17, 17, 17,
                                16, 16, 16, 16, 16))
q
ggsave(plot = q, device = pdf, filename = "results/plots/2DtSNE_Combined_v2.pdf", width = 8, height = 6)


# 3D PCA
library(scatterplot3d)
prData <- prcomp(log2(mat))
pca.data <- prData$rotation
pca.data <- data.frame(pca.data)[1:6]
pca.data <- cbind(pca.data, sampAnnot)  

# Add shape and color
pca.data$shape <- ifelse(pca.data$Type == "RNA-Seq", 16, 17)
pca.data$color[pca.data$Subgroup == "Group3"] <- rgb(red = 248/255, green = 118/255, blue = 109/255, alpha = 0.6)
pca.data$color[pca.data$Subgroup == "Group4"] <- rgb(red = 124/255, green = 174/255, blue = 0/255, alpha = 0.6)
pca.data$color[pca.data$Subgroup == "WNT"] <- rgb(red = 0/255, green = 191/255, blue = 196/255, alpha = 0.6)
pca.data$color[pca.data$Subgroup == "SHH"] <- rgb(red = 199/255, green = 124/255, blue = 255/255, alpha = 0.6)
pca.data$color[pca.data$Subgroup == "Unknown"] <- rgb(red = 0/255, green = 0/255, blue = 0/255, alpha = 0.6)

fname <- paste0('results/plots/3DPCA_Combined.pdf')
pdf(file = fname, width = 12, height = 10)
s3d <- scatterplot3d(pca.data[,"PC1"], pca.data[,"PC2"], pca.data[,"PC3"], 
                     xlab = "PC1", ylab = "PC2", zlab = "PC3", 
                     color = as.character(pca.data[,"color"]), 
                     pch = pca.data[,"shape"], 
                     cex.symbols = 3, main = "PC1-PC3") 
tmpLegend <- unique(pca.data[,c("Type", "shape", "Subgroup", "color")])
tmpLegend[,"Group"] <- paste0(tmpLegend[,"Type"],"_",tmpLegend[,"Subgroup"])
reorder.cols <- match(tmpLegend$Group, c("RNA-Seq_Group3", "RNA-Seq_Group4", "RNA-Seq_SHH", "RNA-Seq_WNT", 
                                         "Microarray_Group3", "Microarray_Group4", "Microarray_SHH", "Microarray_WNT", "Microarray_Unknown"))
tmpLegend <- tmpLegend[reorder.cols,]
legend(x="right", 
       pch = tmpLegend[,"shape"], 
       col = as.character(tmpLegend[,"color"]), 
       legend = tmpLegend[,"Group"], cex = 1, 
       inset = 0.07)
dev.off()

# PC2-4
fname <- paste0('results/plots/3DPCA_Combined_v2.pdf')
pdf(file = fname, width = 12, height = 10)
s3d <- scatterplot3d(pca.data[,"PC2"], pca.data[,"PC3"], pca.data[,"PC4"], 
                     xlab = "PC2", ylab = "PC3", zlab = "PC4", 
                     color = as.character(pca.data[,"color"]), 
                     pch = pca.data[,"shape"], 
                     cex.symbols = 3, main = "PC2-PC4") 
tmpLegend <- unique(pca.data[,c("Type", "shape", "Subgroup", "color")])
tmpLegend[,"Group"] <- paste0(tmpLegend[,"Type"],"_",tmpLegend[,"Subgroup"])
reorder.cols <- match(tmpLegend$Group, c("RNA-Seq_Group3", "RNA-Seq_Group4", "RNA-Seq_SHH", "RNA-Seq_WNT", 
                                         "Microarray_Group3", "Microarray_Group4", "Microarray_SHH", "Microarray_WNT", "Microarray_Unknown"))
tmpLegend <- tmpLegend[reorder.cols,]
legend(x="right", 
       pch = tmpLegend[,"shape"], 
       col = as.character(tmpLegend[,"color"]), 
       legend = tmpLegend[,"Group"], cex = 1, 
       inset = 0.07)
dev.off()

# PC 3-5
fname <- paste0('results/plots/3DPCA_Combined_v3.pdf')
pdf(file = fname, width = 12, height = 10)
s3d <- scatterplot3d(pca.data[,"PC3"], pca.data[,"PC4"], pca.data[,"PC5"], 
                     xlab = "PC3", ylab = "PC4", zlab = "PC5", 
                     color = as.character(pca.data[,"color"]), 
                     pch = pca.data[,"shape"], 
                     cex.symbols = 3, main = "PC3-PC5") 
tmpLegend <- unique(pca.data[,c("Type", "shape", "Subgroup", "color")])
tmpLegend[,"Group"] <- paste0(tmpLegend[,"Type"],"_",tmpLegend[,"Subgroup"])
reorder.cols <- match(tmpLegend$Group, c("RNA-Seq_Group3", "RNA-Seq_Group4", "RNA-Seq_SHH", "RNA-Seq_WNT", 
                                         "Microarray_Group3", "Microarray_Group4", "Microarray_SHH", "Microarray_WNT", "Microarray_Unknown"))
tmpLegend <- tmpLegend[reorder.cols,]
legend(x="right", 
       pch = tmpLegend[,"shape"], 
       col = as.character(tmpLegend[,"color"]), 
       legend = tmpLegend[,"Group"], cex = 1, 
       inset = 0.07)
dev.off()
# tmpLegendShape <- unique(pca.data[,c("Type", "shape")])
# tmpLegendColor <- unique(pca.data[,c("Subgroup", "color")])
# legend(x = 1, y = 4, 
#        pch = tmpLegendShape[,"shape"], 
#        legend = tmpLegendShape[,"Type"], cex = 0.7, xpd = TRUE, inset = -0.3)
# legend(x = 1, y = 2,
#        fill = tmpLegendColor[,"color"], 
#        legend = tmpLegendColor[,"Subgroup"], cex = 0.6,  xpd = TRUE, inset = -0.2)

