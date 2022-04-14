library(tidyverse)
library(preprocessCore)
library(survival)
library(survminer)
library(DESeq2)
library(ComplexHeatmap)

set.seed(2020)

setwd("/Users/ferrenaa/Documents/deyou/target/")




#tpm or counts
tpm = T


if(tpm == T){
  l <- readRDS('data/parsed/parsed_tpm.rds')
} else{
  l <- readRDS('data/parsed/parsed.rds')
}

#read the list:
# names(l)
# 1           gem
# 2 samp_case_ids
# 3 clinical_data
# 4      metadata
# 5  ensembl_hugo


gem  <- l[[1]]
sdf  <- l[[2]]
clin <- l[[3]]
sampmat <- l[[4]]
genes <- l[[5]]

rm(l, sampmat, genes)


#remove missing survival
clin <- clin[!is.na(clin$`Event Free Survival Time in Days`),]
sdf <- sdf[sdf$cases %in% clin$`TARGET USI`,]
gem <- gem[,colnames(gem) %in% sdf$samps]


#check lib size
libsize <- data.frame(samp = colnames(gem),
                      libsize = colSums(gem))



ggplot(libsize, aes(x = libsize))+
  geom_histogram(col = 'black', fill='steelblue')+
  scale_x_log10()



#remove extreme low TPM samples
badsamps <- c('TARGET-40-PALFYN-01A-01R',
              'TARGET-40-PASKZZ-01A-01R')


sdf <- sdf[!(sdf$samps %in% badsamps),]
gem <- gem[,colnames(gem) %in% sdf$samps]
clin <- clin[clin$`TARGET USI` %in% sdf$cases,]


#norm: quant norm for TPM, DESeq2 size factors if counts

if(tpm == T){
  
  #quantile norm; 
  gemt <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(gem)))
  rownames(gemt) <- rownames(gem)
  colnames(gemt) <- colnames(gem)
  
  gem <- gemt
  rm(gemt)
  
} else{
  
  gem <- t( t(gem) / DESeq2::estimateSizeFactorsForMatrix(gem) )
  
}






#try for DKOAA module score
res <- read.csv('~/Documents/deyou/dko_dkoaa_october2020/start_to_end/results_figures/DEG_filtered_padj0.05.csv')[,1:8]
res <- res[res$padj < 0.05,]

#get mouse homologs

library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

convertMouseGeneList <- function(x){
  require("biomaRt")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), 
                   filters = "mgi_symbol", 
                   values = x , 
                   mart = mouse, 
                   attributesL = c("hgnc_symbol"), 
                   martL = human, 
                   uniqueRows=T)
  
  genesV2 
}

genes <- res$Gene_name
hugo <- convertMouseGeneList(genes)

#remove duplicates...
hugo <- hugo[!duplicated(hugo[,1]),]
hugo <- hugo[!duplicated(hugo[,2]),]
res <- res[res$Gene_name %in% hugo$MGI.symbol,]

#match order
hugo <- hugo[match(res$Gene_name, hugo$MGI.symbol),]

res$hugo <- hugo$HGNC.symbol


dkoaa_overexp  <- res[res$log2FoldChange>0,'hugo']
dkoaa_underexp <- res[res$log2FoldChange<0,'hugo']







#try to scale gem...
gem <- t(scale(t(gem)))





dkoaa_overexp <- dkoaa_overexp[dkoaa_overexp %in% rownames(gem)]



#test module score, look at it

#real module score
# should correlate wll with real genes
# should correlate poorly with fake genes


mod <- scale(modulescore(gem, dkoaa_overexp))


#fake module score; select random genes
# should correlate poor;ly with real genes
# should correlate better with fake genes

fakegenes <- sample(rownames(gem), size = length(dkoaa_overexp), replace = F)
fakemod <- scale(modulescore(gem, fakegenes) )




#first test if module score correlates with real genes
subgem <- gem[rownames(gem) %in% dkoaa_overexp,]


#correlate module score with each gene
corres <- apply(subgem, MARGIN = 1, FUN = function(x){
  cor(x, mod)
})


hist(corres, xlim=c(-1,1))




## next test if module score correlates with fake genes
fakesubgem <- gem[rownames(gem) %in% fakegenes,]


#correlate module score with each gene
fakecorres <- apply(fakesubgem, MARGIN = 1, FUN = function(x){
  cor(x, mod)
})
hist(fakecorres,  xlim=c(-1,1))


#compare them


plot(corres, fakecorres, xlim=c(-1,1), ylim=c(-1,1))
abline(a = 0, b = 1, col='red')

t.test(corres, fakecorres)




#fake module score: should correlate poorly with real genes
corres_fakemod <- apply(subgem, MARGIN = 1, FUN = function(x){
  cor(x, fakemod)
})
hist(corres_fakemod,  xlim=c(-1,1))



#fake module score: should correlate poorly with real genes
fakecorres_fakemod <- apply(fakesubgem, MARGIN = 1, FUN = function(x){
  cor(x, fakemod)
})
hist(fakecorres_fakemod,  xlim=c(-1,1))



plot(corres_fakemod, fakecorres_fakemod, xlim=c(-1,1), ylim=c(-1,1))
abline(a = 0, b = 1, col='red')

t.test(corres_fakemod, fakecorres_fakemod)





mean(corres)





#plots
par(mfrow=c(2,1))
hist(corres, xlim=c(-1,1))
hist(fakecorres,  xlim=c(-1,1))

hist(corres_fakemod,  xlim=c(-1,1))
hist(fakecorres_fakemod,  xlim=c(-1,1))
dev.off()

plot(corres, fakecorres, xlim=c(-1,1), ylim=c(-1,1))
abline(a = 0, b = 1, col='red')
t.test(corres, fakecorres)


plot(corres_fakemod, fakecorres_fakemod, xlim=c(-1,1), ylim=c(-1,1))
abline(a = 0, b = 1, col='red')
t.test(corres_fakemod, fakecorres_fakemod)






#randomly sample
ntests <- 









#correlation: all genes + module score
subgem <- rbind(mod, gem[rownames(gem) %in% dkoaa_overexp,])
rownames(subgem)[1] <- 'modulescore'


corres <- cor(t(subgem))


library(ComplexHeatmap)

labs <- c("ModuleScore", rep('', nrow(subgem) - 1))
Heatmap(corres, cluster_rows = F, cluster_columns = F, 
        row_labels = labs, column_labels = labs)


#cor for fake genes
fakegens <- sample(rownames(gem), size = length(dkoaa_overexp), replace = F)
subgem <- gem[rownames(gem) %in% fakegenes,]
rownames(subgem)[1] <- 'modulescore'


corres <- cor(t(subgem))


library(ComplexHeatmap)

labs <- c("ModuleScore", rep('', nrow(subgem) - 1))
Heatmap(corres, cluster_rows = F, cluster_columns = F, 
        row_labels = labs, column_labels = labs)







#correlation: all genes + module score
fakegenes <- sample(rownames(gem), size = length(dkoaa_overexp), replace = F)
fakemod <- modulescore(gem, fakegenes)


subgem <- rbind(fakemod, gem[rownames(gem) %in% dkoaa_overexp,])
rownames(subgem)[1] <- 'modulescore'


corres <- cor(t(subgem))


library(ComplexHeatmap)

labs <- c("ModuleScore", rep('', nrow(subgem) - 1))
Heatmap(corres, cluster_rows = F, cluster_columns = F, 
        row_labels = labs, column_labels = labs)


#cor for fake genes

subgem <- gem[rownames(gem) %in% fakegenes,]
rownames(subgem)[1] <- 'modulescore'


corres <- cor(t(subgem))


library(ComplexHeatmap)

labs <- c("ModuleScore", rep('', nrow(subgem) - 1))
Heatmap(corres, cluster_rows = F, cluster_columns = F, 
        row_labels = labs, column_labels = labs)






