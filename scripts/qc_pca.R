library(tidyverse)
library(preprocessCore)
library(survival)
library(survminer)
library(DESeq2)

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

rm(gemt, gem_hvgs, rowvars, gem_scale)

plot(pca$sdev)

#get embeddings
emb <- as.data.frame(pca$x)

#get genes
pcagenes <- pca$rotation
head(sort(pcagenes[,1],decreasing = T))
head(sort(pcagenes[,2],decreasing = T))

libsize <- libsize[libsize$samp %in% sdf$samps,]
emb$libsize <- libsize$libsize

#plot PCs
ggplot(emb, aes(PC1, PC2, col=libsize))+
  geom_point()+
  theme_light()




#color by some gene
coolgenes <- c('HOXA9', 'COX7C', 'PTPRZ1')


coolgenesgem <- gem[rownames(gem) %in% coolgenes,]


pl <- list()
for(gene in coolgenes){
  
  #emb$gene <- scale(as.numeric(coolgenesgem[gene,]))
  emb$gene <- as.numeric(coolgenesgem[gene,])
  
  g <- ggplot(emb, aes(PC1, PC2, col=gene))+
    geom_point()+
    theme_light()+
    scale_color_distiller(palette = 'RdBu') + labs(color = gene)
  
  pl[[gene]] <- g
}

cowplot::plot_grid(plotlist = pl)

  
### do that later ###



gene <- 'SKP2'
genevec <- as.numeric(gem[gene,])
names(genevec) <- colnames(gem)

df <- data.frame(samps = sdf$samps,
                 cases = sdf$cases,
                 gene = genevec,
                 gene_dich = ifelse(genevec >= median(genevec), yes = 'hi', no = 'lo'), 
                 surv = clin$`Event Free Survival Time in Days` / 30,
                 status = clin$`Vital Status`)

df$gene_dich <- factor(df$gene_dich, levels = c('lo', 'hi'))
ggplot(df, aes(group = gene_dich, y = log2(gene), x = gene_dich))+
  geom_violin()+geom_jitter(width = 0.1)


#recode: 1 = censored, 2 = dead
df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)








df <- data.frame(samps = sdf$samps,
                 cases = sdf$cases,
                 gene = genevec,
                 gene_dich = ifelse(genevec >= median(genevec), yes = 'hi', no = 'lo'), 
                 surv = clin$`Overall Survival Time in Days`,
                 status = clin$`Vital Status`)


df$surv <- df$surv / 30



df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)

fit <- survfit(Surv(surv, status_code) ~ 1, data = df)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           linetype = "strata", 
           surv.median.line = "hv", 
           ggtheme = theme_light())+ xlab('Time (Months)')






df <- data.frame(samps = sdf$samps,
                 cases = sdf$cases,
                 gene = genevec,
                 gene_dich = ifelse(genevec >= median(genevec), yes = 'hi', no = 'lo'), 
                 surv = clin$`Overall Survival Time in Days`,
                 status = clin$`Vital Status`)


df$surv <- df$surv / 365



df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)

fit <- survfit(Surv(surv, status_code) ~ 1, data = df)

ggsurvplot(fit,
           pval = TRUE, conf.int = F,
           linetype = "strata", 
           #surv.median.line = "hv", 
           #risk.table = T, cumevents = T, cumcensor = T, tables.height = 0.1,
           ggtheme = theme_light())+ xlab('Time (Years)')












### HOXA9 and other genes ###

gene <- 'HOXA9'
gene <- 'SKP2'
gene <- 'MMP9'
gene <- 'CKS1B'
gene <- 'CDKN1B'
gene <- 'TP53'
gene <- 'RB1'
gene <- 'RHOXF1'
gene <- 'KRT7'

### overall surv ###
time <- "Overall Survival Time in Days"

genevec <- as.numeric(gem[gene,])
names(genevec) <- colnames(gem)

df <- data.frame(samps = sdf$samps,
                 cases = sdf$cases,
                 gene = genevec,
                 gene_dich = ifelse(genevec >= median(genevec), yes = 'hi', no = 'lo'), 
                 surv = clin[,time],
                 status = clin$`Vital Status`)

df$surv <- df$surv / 30

df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)

fit <- survfit(Surv(surv, status_code) ~ gene_dich, data = df)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           linetype = "strata", 
           title=gene, subtitle=time,
           surv.median.line = "hv", 
           ggtheme = theme_light(), 
           palette = RColorBrewer::brewer.pal('Set1',n=3)[1:2])+ xlab('Time (Months)')






### event free surv ###
time <- "Event Free Survival Time in Days" 
gene <- 'SKP2'
genevec <- as.numeric(gem[gene,])
names(genevec) <- colnames(gem)

df <- data.frame(samps = sdf$samps,
                 cases = sdf$cases,
                 gene = genevec,
                 gene_dich = ifelse(genevec >= median(genevec), yes = 'hi', no = 'lo'), 
                 surv = clin[,time],
                 status = clin$`Vital Status`)

df$surv <- df$surv / 30

df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)

fit <- survfit(Surv(surv, status_code) ~ gene_dich, data = df)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           linetype = "strata", 
           title=gene, subtitle=time,
           surv.median.line = "hv", 
           ggtheme = theme_light(), 
           palette = RColorBrewer::brewer.pal('Set1',n=3)[1:2])+ xlab('Time (Months)')



### "Time to first relapse in days"
time <- "Time to first relapse in days" 
gene <- 'SKP2'
genevec <- as.numeric(gem[gene,])
names(genevec) <- colnames(gem)

df <- data.frame(samps = sdf$samps,
                 cases = sdf$cases,
                 gene = genevec,
                 gene_dich = ifelse(genevec >= median(genevec), yes = 'hi', no = 'lo'), 
                 surv = clin[,time],
                 status = clin$`Vital Status`)

df$surv <- df$surv / 30

df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)

fit <- survfit(Surv(surv, status_code) ~ gene_dich, data = df)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           linetype = "strata", 
           title=gene, subtitle=time,
           surv.median.line = "hv", 
           ggtheme = theme_light(), 
           palette = RColorBrewer::brewer.pal('Set1',n=3)[1:2])+ xlab('Time (Months)')







### "Time to death in days" 
time <- "Time to death in days"  
gene <- 'SKP2'
genevec <- as.numeric(gem[gene,])
names(genevec) <- colnames(gem)

df <- data.frame(samps = sdf$samps,
                 cases = sdf$cases,
                 gene = genevec,
                 gene_dich = ifelse(genevec >= median(genevec), yes = 'hi', no = 'lo'), 
                 surv = clin[,time],
                 status = clin$`Vital Status`)

df$surv <- df$surv / 30

df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)

fit <- survfit(Surv(surv, status_code) ~ gene_dich, data = df)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           linetype = "strata", 
           title=gene, subtitle=time,
           surv.median.line = "hv", 
           ggtheme = theme_light(), 
           palette = RColorBrewer::brewer.pal('Set1',n=3)[1:2])+ xlab('Time (Months)')






#loop thru all genes...

genes <- rownames(gem)

#progress bar
total = length(genes) - 1
pb <- txtProgressBar(min = 0, max = total, style = 3)




reslist <- list()
for(genedex in c(1:length(genes)) ){
  
  gene <- genes[genedex]
  
  
  genevec <- as.numeric(gem[gene,])
  names(genevec) <- colnames(gem)
  
  
  outdf <- tryCatch({
    df <- data.frame(samps = sdf$samps,
                     cases = sdf$cases,
                     gene = genevec,
                     gene_dich = ifelse(genevec >= median(genevec), yes = 1, no = 0), 
                     surv = clin$`Overall Survival Time in Days`,
                     status = clin$`Vital Status`)
    
    #recode: 1 = censored, 2 = dead
    df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)
    
    chisq <- survdiff(Surv(surv, status_code) ~ gene_dich, data = df)$chisq
    
    p <- pchisq(chisq, df=1, lower.tail=FALSE)
    
    outdf <- data.frame(gene =gene, p = p)
    outdf
  }, error = function(e){
    outdf <- data.frame(gene = gene, p = NA)
    outdf
  })
  
  reslist[[genedex]] <- outdf
  setTxtProgressBar(pb, genedex)
  
}

beepr::beep()



res <- bind_rows(reslist)
res <- res[!is.na(res$p),]
res$fdr <- p.adjust(res$p)


head(res[order(res$fdr),])
head(res[order(res$p),])


gene <- 'HAVCR2'
gene <- 'SIRPA'
gene <- 'GATA3'
gene <- 'IRF2BPL'
gene <- 'MMD'
gene <- 'MAFK'

gene <- 'SIRPA'
gene <- 'IRF2BPL'
gene <- 'MT1A'



### overall surv ###
time <- "Overall Survival Time in Days"
gene <- 'HAVCR2'
genevec <- as.numeric(gem[gene,])
names(genevec) <- colnames(gem)

df <- data.frame(samps = sdf$samps,
                 cases = sdf$cases,
                 gene = genevec,
                 gene_dich = ifelse(genevec >= median(genevec), yes = 'hi', no = 'lo'), 
                 surv = clin[,time],
                 status = clin$`Vital Status`)

df$surv <- df$surv / 30

df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)

fit <- survfit(Surv(surv, status_code) ~ gene_dich, data = df)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           linetype = "strata", 
           title=gene, subtitle=time,
           surv.median.line = "hv", 
           ggtheme = theme_light(), 
           palette = RColorBrewer::brewer.pal('Set1',n=3)[1:2])+ xlab('Time (Months)')




### overall surv ###
time <- "Overall Survival Time in Days"
gene <- 'MT1A'
genevec <- as.numeric(gem[gene,])
names(genevec) <- colnames(gem)

df <- data.frame(samps = sdf$samps,
                 cases = sdf$cases,
                 gene = genevec,
                 gene_dich = ifelse(genevec >= median(genevec), yes = 'hi', no = 'lo'), 
                 surv = clin[,time],
                 status = clin$`Vital Status`)

df$surv <- df$surv / 30

df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)

fit <- survfit(Surv(surv, status_code) ~ gene_dich, data = df)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           linetype = "strata", 
           title=gene, subtitle=time,
           surv.median.line = "hv", 
           ggtheme = theme_light(), 
           palette = RColorBrewer::brewer.pal('Set1',n=3)[1:2])+ xlab('Time (Months)')













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


set.seed(2020)


### overexp ###

title <- 'Genes Overepxressed in DKOAA'
mod <- modulescore(gem, dkoaa_overexp)
time <- "Overall Survival Time in Days" 

df <- data.frame(samps = sdf$samps,
                 cases = sdf$cases,
                 mod = mod,
                 mod_dich = ifelse(mod >= median(mod), yes = 'hi', no = 'lo'), 
                 surv = clin[,time],
                 status = clin$`Vital Status`)

df$surv <- df$surv / 30

df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)

fit <- survfit(Surv(surv, status_code) ~ mod_dich, data = df)

ggsurvplot(fit,
           pval = TRUE, conf.int = F,
           linetype = "strata", 
           title=title, subtitle=time,
           #surv.median.line = "hv", 
           ggtheme = theme_light(), 
           palette = c('Black','Grey50')
           ) + 
  xlab('Overall survival (Months)')




### underexp ###
title <- 'Genes underepxressed in DKOAA'
mod <- modulescore(gem, dkoaa_underexp)
time <- "Overall Survival Time in Days" 

df <- data.frame(samps = sdf$samps,
                 cases = sdf$cases,
                 mod = mod,
                 mod_dich = ifelse(mod >= median(mod), yes = 'hi', no = 'lo'), 
                 surv = clin[,time],
                 status = clin$`Vital Status`)

df$surv <- df$surv / 30

df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)

fit <- survfit(Surv(surv, status_code) ~ mod_dich, data = df)

ggsurvplot(fit,
           pval = TRUE, conf.int = F,
           linetype = "strata", 
           title=title, subtitle=time,
           #surv.median.line = "hv", 
           ggtheme = theme_light(), 
           palette = c('Black','Grey50')
) + 
  xlab('Overall survival (Months)')

  





