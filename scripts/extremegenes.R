library(tidyverse)
library(preprocessCore)
library(survival)
library(survminer)
library(DESeq2)

set.seed(2020)

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



### read in cibersort ###
tm <- read.csv('data/cibersortresults/CIBERSORTx_Job2_Results.csv')
ith <- read.csv('data/cibersortresults/CIBERSORTx_Job6_Results.csv')

bc16 <- read.csv('data/cibersortresults/bc16/CIBERSORTx_Job12_Results.csv')

#add it in
tm <- tm[match(sdf$samps, tm$Mixture),]
tm$Microenv <- rowSums(tm[,3:10]) #total microenv, sum non-tumor
ith <- ith[match(sdf$samps, ith$Mixture),]

bc16 <- bc16[match(sdf$samps, bc16$Mixture),]
bc16$Microenv <- rowSums(bc16[,3:7])

#exclude the uneccessary columns
excludecols <- c('Mixture', 'P.value', 'Correlation', 'RMSE', 'Absolute.score..sig.score.')
tmgood <- tm[,!(colnames(tm) %in% excludecols)]
ithgood <- ith[,!(colnames(ith) %in% excludecols)]
bc16good <- bc16[,!(colnames(bc16) %in% excludecols)]

#for bc16, append bc16 to colnames
colnames(bc16good) <- paste0('BC16_', colnames(bc16good))

#add to clin, but get clin names first
origclin <- colnames(clin)

clin <- cbind(clin,
              tmgood/tm$Absolute.score..sig.score.,
              ithgood/ith$Absolute.score..sig.score.,
              bc16good/bc16$Absolute.score..sig.score.)

### see the correlation of human and mouse cibersort results ###


shared <- data.frame(mouse_tumor = clin$Tumor, human_tumor = clin$BC16_Tumor,
                     mouse_macrophage = clin$Monocyte.Macrophage, human_macrophage = clin$BC16_Macrophage,
                     mouse_osteoclast = clin$Osteoclast, human_osteoclast = clin$BC16_Osteoclast,
                     mouse_tcell = clin$Tcell, human_tcell = clin$BC16_Tcell,
                     mouse_endothelial = clin$Endothelial, human_endothelial = clin$BC16_Endothelial,
                     mouse_fibroblast = clin$MSC.Fibroblast, human_fibroblast = clin$BC16_Fibroblast,
                     mouse_microenv = clin$Microenv, human_microenv = clin$BC16_Microenv
)


sharedcelltypes <- unique(substring(colnames(shared), first = 7))

plotlist <- list()
for(ct in sharedcelltypes){
  ctdf <- shared[,grep(ct, colnames(shared))]
  
  exactlabels <- colnames(ctdf)
  colnames(ctdf) <- c('mouse', 'human')
  
  cres <- cor.test(ctdf[,1], ctdf[,2])
  nicep <- ifelse(test = cres$p.value < 0.01, formatC(cres$p.value, format='e', digits=2), 
                  no = round(cres$p.value, 2))
  
  
  
  g <- ggplot(ctdf, aes(mouse, human))+
    geom_point(size=3)+
    labs(title = paste0(exactlabels[1], ' vs ', exactlabels[2]),
         subtitle = paste0('Pearson cor r = ', round(cres$estimate, 2),
                           '\nP (cor.test) = ', nicep))+
    ylab(exactlabels[1])+ xlab(exactlabels[2])+
    theme_linedraw() 
  
  plotlist[[ct]] <- g
  
}

patchwork::wrap_plots(plotlist)

patchwork::wrap_plots(plotlist[2:6])

plotlist[[1]] / plotlist[[7]]




#plot the scores first
means <- colMeans(tmgood)
#exclude microenv...
means <- means[names(means) != 'Microenv']
means <- means / sum(means)
meandf <- data.frame(names = names(means), mean=means)

meandf$names <- factor(meandf$names, levels=names( sort(means)))


set.seed(420)
wholepal <- sample( scales::hue_pal()(nrow(meandf)), replace = F)
tmprop <- ggplot(meandf, aes(x='',y=mean, fill=names))+
  geom_bar(stat='identity', position = 'stack')+
  scale_fill_manual(values =  wholepal)


#plot the scores first
means <- colMeans(ithgood) / sum(colMeans(ithgood))
meandf <- data.frame(names = names(means), mean=means)

meandf$names <- factor(meandf$names, levels=names( sort(means)))


set.seed(69)
tumorpal <- sample( scales::hue_pal()(nrow(meandf)), replace = F)
ithprop <- ggplot(meandf, aes(x='',y=mean, fill=names))+
  geom_bar(stat='identity', position = 'stack')+
  scale_fill_manual(values =  tumorpal)




#plot the scores first
means <- colMeans(bc16good)
#exclude microenv...
means <- means[names(means) != 'BC16_Microenv']
means <- means / sum(means)
meandf <- data.frame(names = names(means), mean=means)

meandf$names <- factor(meandf$names, levels=names( sort(means)))


set.seed(420)
wholepal <- sample( scales::hue_pal()(nrow(meandf)), replace = F)
bcprop <- ggplot(meandf, aes(x='',y=mean, fill=names))+
  geom_bar(stat='identity', position = 'stack')+
  scale_fill_manual(values =  wholepal)






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



#check lib size
libsize_post <- data.frame(samp = colnames(gem),
                           libsize = colSums(gem))



ggplot(libsize_post, aes(x = libsize))+
  geom_histogram(col = 'black', fill='steelblue')



#PCA

numhvgs <- 2000

#transform; sqrt or log. can also try VST or rlog...
gemt <- sqrt(gem)


#get HVGs
rowvars <- apply(gemt, 1, var)

#select top hvgs
rowvars <- sort(rowvars, decreasing = T)[1:numhvgs]

#get a small mat
gem_hvgs <- gemt[rownames(gemt) %in% names(rowvars),]

#scale small mat
gem_scale <- t(scale(t(gem_hvgs)))




#pca
pca <- prcomp(t(gem_scale))

#rm(gemt, gem_hvgs, rowvars, gem_scale)

plot(pca$sdev)

#get embeddings
emb <- as.data.frame(pca$x)

#get genes
pcagenes <- pca$rotation
head(sort(pcagenes[,1],decreasing = T))
head(sort(pcagenes[,2],decreasing = T))

libsize <- libsize[libsize$samp %in% sdf$samps,]
emb$libsize <- libsize$libsize

emb$Tumor <- clin$Tumor
emb$Microenv <- clin$Microenv

emb$BC16_Tumor <- clin$BC16_Tumor
emb$BC16_Microenv <- clin$BC16_Microenv

#plot PCs
ggplot(emb, aes(PC1, PC2, col=libsize))+
  geom_point()+
  theme_light()+
  scale_color_distiller(palette = 'Purples')









### check extreme genes ###
numhvgs <- 500

#transform; sqrt or log. can also try VST or rlog...
gemt <- t(apply(gem, 1, function(x){log2(x+1)}))

#get HVGs
rowvars <- apply(gemt, 1, var)

#select top hvgs
rowvarstop <- sort(rowvars, decreasing = T)[1:numhvgs]

#make mean0var plot
rowmeans <- rowMeans(gemt)
rowvars <- rowvars

gdf <- data.frame(rowmeans, rowvars)
  
gdf$hvg <- 'Non-HVG'
gdf[rownames(gdf) %in% names(rowvarstop), "hvg"] <- "HVG"

top <- head(gdf[order(gdf$rowvars, decreasing = T),], n=50)
top$repel <- rownames(top)


ggplot(gdf, aes(x=rowmeans, y=rowvars, col=hvg))+
  geom_point(size = 0.5)+
  ggrepel::geom_text_repel(data = top, inherit.aes = F,aes(x=rowmeans, y=rowvars,label=repel), size= 3)+
  scale_color_brewer(palette = 'Set1')+
  theme_light()+
  labs(x='Mean Log2 Expression', y = 'Variance of Log2 Gene Expression')




#get a small mat
gem_hvgs <- gemt[rownames(gemt) %in% names(rowvarstop),]


#scale small mat
gem_scale <- t(scale(t(gem_hvgs)))

tmpgem <- gem_hvgs

#heatmap
colors <- rev(colorRampPalette( brewer.pal(9, "RdYlBu"))(255))
colors <- circlize::colorRamp2(seq(min(tmpgem), max(tmpgem), length = length(colors)), colors)


Heatmap(tmpgem, name = 'HVGs',
        col = colors, 
        row_labels = rep("", nrow(gem_hvgs)),
        column_labels = rep("", ncol(gem_hvgs))
        )



