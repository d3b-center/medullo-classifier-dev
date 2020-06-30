# Author: Komal S. Rathi
# Date: 10/19/2019
# Function:  
# Accuracy comparison of MB classifier with MM2S package
# SuppFigure3

library(MM2S)
library(MM2Sdata)
library(preprocessCore)
library(medulloPackage)
library(plyr)
library(xlsx)

source('R/utils/pubTheme.R')

# MM2S package uses ENTREZ gene ids whereas we use gene symbols

# Training set
load('data/loadedGSE_37418.RData')
geneAnnot <- read.delim("data/GPL570-55999.txt", skip=16)
geneAnnot <- geneAnnot[,c("ID", "ENTREZ_GENE_ID")]
exprs_37418_Tmp <- normalize.quantiles(as.matrix(exprs_37418))
rownames(exprs_37418_Tmp) <- rownames(exprs_37418)
colnames(exprs_37418_Tmp )<- colnames(exprs_37418)
exprs_37418 <- exprs_37418_Tmp
exprs_37418 <- data.frame(exprs_37418)
exprs_37418[,"Max"] <- apply(exprs_37418, FUN=max, MARGIN=1)
exprs_37418[,"Probe"] <- rownames(exprs_37418)
exprs_37418 <- merge(exprs_37418, geneAnnot, by.x="Probe", by.y="ID")
exprs_37418 <- exprs_37418[order(-exprs_37418[,"Max"]),]

# Remove rows without gene symbols & make sure one symbol per row
exprs_37418 <- exprs_37418[exprs_37418[,"ENTREZ_GENE_ID"]!="",]
exprs_37418 <- exprs_37418[!duplicated(exprs_37418[,"ENTREZ_GENE_ID"]),]
rownames(exprs_37418) <- exprs_37418[,"ENTREZ_GENE_ID"]
exprs_37418 <- exprs_37418[-1]
exprs_37418 <- exprs_37418[1:(ncol(exprs_37418)-2)]

# observed subtypes
observed_37418 <- annot_37418[,c('geo_accession','characteristics_ch1')]
observed_37418$characteristics_ch1 <- gsub('subgroup: ', '', observed_37418$characteristics_ch1)
observed_37418$characteristics_ch1 <- gsub("G3", "Group3", observed_37418$characteristics_ch1)
observed_37418$characteristics_ch1 <- gsub("G4", "Group4", observed_37418$characteristics_ch1)
observed_37418$characteristics_ch1 <- gsub("SHH OUTLIER", "SHH", observed_37418$characteristics_ch1)
observed_37418$characteristics_ch1 <- gsub("U", "Unknown", observed_37418$characteristics_ch1)

# Predict subtypes using MM2S package
pred_37418 <- MM2S.human(InputMatrix = exprs_37418,
                         parallelize = 4, 
                         seed = 12345, tempdir())
pred_37418 <- as.data.frame(pred_37418$MM2S_Subtype)

# merge observed and predicted
res <- observed_37418 %>%
  inner_join(pred_37418, by = c("geo_accession" = "SampleName"))

# calculate accuracy
res <- calcStats(myClassActual = res$characteristics_ch1,  myClassPred = res$MM2S_Prediction)
accuracy <- trimws(res[[2]]$stats[[1]])

# now perform accuracy test on all 15 test datasets
load('data/expr.RData')
load('data/meta.RData')
colnames(actual)[3] <- 'observed'
study.ct <- plyr::count(actual$study)

# convert gene symbols to Entrez ids
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values = rownames(dat), bmHeader = T, mart = mart)
colnames(mapping) = c("hgnc_symbol","ncbi_gene_id")
mapping <- mapping %>%
  filter(!is.na(ncbi_gene_id))

exprTest <- mapping %>%
  inner_join(dat %>%
               rownames_to_column('hgnc_symbol'), by = 'hgnc_symbol') %>%
  dplyr::select(-c(hgnc_symbol)) %>%
  mutate(means = rowMeans(dplyr::select(.,-ncbi_gene_id))) %>%
  arrange(desc(means)) %>% # arrange decreasing by means
  distinct(ncbi_gene_id, .keep_all = TRUE) %>% # keep the ones with greatest mean value. If ties occur, keep the first occurencce
  dplyr::select(-c(means)) %>%
  unique() %>%
  remove_rownames() %>%
  column_to_rownames('ncbi_gene_id')

#x <- actual[which(actual$study ==  "EMTAB292"),]
#dat <- exprTest
run.mm2s <- function(meta, dat, type  = c("pred", "stats")){
  # first predict
  print(paste0('Study: ', unique(meta$study)))
  rownames(meta) <- meta$id
  dat.sub <- dat[,rownames(meta)]
  
  # predict using MM2S
  pred <- MM2S.human(InputMatrix = dat.sub,
                     parallelize = 4, 
                     seed = 12345, tempdir())
  pred <- as.data.frame(pred$MM2S_Subtype)
  total <- meta %>%
    inner_join(pred, by = c("id" = "SampleName"))
  if(type == "pred") {
    return(total)
  }
  
  # test test
  res <- calcStats(myClassActual = total$observed,  myClassPred = total$MM2S_Prediction)
  accuracy <- res[[2]]$stats[[1]]
  sens_spec <- res[[3]][1:2]
  rownames(sens_spec) <- gsub("Class: ","", rownames(sens_spec))
  sens_spec <- melt(as.matrix(sens_spec))
  df <- data.frame(var = paste0(sens_spec$Var1,"_", sens_spec$Var2), val = sens_spec$value)
  df <- rbind(df, data.frame(var = "Accuracy", val = accuracy))
  return(df)
}

# save stats
stats <- ddply(actual, .variables = "study", .fun = function(x) run.mm2s(meta = x, dat, type = "stats"))
big.res <- dcast(stats, study~var, value.var = 'val')
save(big.res, study.ct, file = 'results/MM2S_GSE124814_split_results.RData')

# save predictions
pred <- ddply(actual, .variables = "study", .fun = function(x) run.mm2s(meta = x, dat, type = "pred"))
save(pred, file = 'results/MM2S_GSE124814_split_results_predictions.RData')
studies <- unique(pred$study)
for(i in 1:length(studies)){
  tmp <- pred[which(pred$study == studies[i]),]
  write.xlsx(tmp, file = "results/tables/MM2S_GSE124814_predictions.xlsx", sheetName = studies[i], row.names = F, append = TRUE)
}

# output accuracy metrics in Table3
res <- merge(study.ct, big.res, by.x = 'x', by.y = 'study')
colnames(res)[1:2] <- c('Study','Sample_Size')
res[res == "NA%"] <- NA
write.xlsx(x = res, file = 'results/tables/MM2S_GSE124814_metrics.xlsx', row.names = F)

format.res <- function(big.res){
  big.res <- big.res[which(big.res$Accuracy != "NaN%"),]
  accuracy <- big.res[,c('study','Accuracy')]
  cor <- as.numeric(gsub('%','',accuracy$Accuracy))
  studies <- accuracy$study
  n <- study.ct[study.ct$x %in% studies,'freq']
  
  # format to add labels
  for.meta <- data.frame(studies, n, accuracy = cor)
  rownames(for.meta) <- for.meta$studies
  for.meta$label = paste0(for.meta$studies, '\n(n = ',for.meta$n,')')
  for.meta <- for.meta[order(for.meta$accuracy, for.meta$n),]
  return(for.meta)
}

# format data
MM2S.res <- format.res(big.res = big.res)
MM2S.res$classifier <- 'MM2S'

# now do a comparison of accuracy, sensitivity and specificity
medulloPackage.res <- get(load('results/GSE124814_split_results.RData'))
medulloPackage.res <- format.res(medulloPackage.res)
medulloPackage.res$classifier <- 'MB Classifier'

# create dodge barplot
plot.dat <- rbind(MM2S.res, medulloPackage.res)
plot.dat$label <- factor(plot.dat$label, levels = medulloPackage.res$label)
plot.dat$accuracy <-  round(plot.dat$accuracy, digits = 1)
p <- ggplot(plot.dat, aes(x = label, y = accuracy, fill = classifier)) +
  geom_bar(position="dodge", stat="identity")  +
  geom_text(aes(label = paste0(accuracy,"%"), vjust = -0.5), size = 2.5, color = "black", position = position_dodge(width = 1)) +
  theme_Publication(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Accuracy (%)") + xlab("") +
  geom_abline(slope = 0, intercept = median(sort(medulloPackage.res$accuracy)),  col = "steelblue", lty = 2, lwd = 0.3) +
  geom_abline(slope = 0, intercept = median(sort(MM2S.res$accuracy)),  col = "lightblue", lty = 2, lwd = 0.3) +
  scale_fill_brewer(palette = "Paired", direction = -1) +
  ggtitle("Accuracy comparison between MB Classifier and MM2S\n15 microarray test datasets")
p
ggsave(p, device = "pdf", filename = 'results/plots/SuppFigure3.pdf', width = 10, height = 6)

