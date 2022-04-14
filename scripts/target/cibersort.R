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

pcam <- ggplot(emb, aes(PC1, PC2, col=Tumor))+
  geom_point(size=3)+
  theme_light()+
  scale_color_distiller(palette='Reds',direction = 1)

pcat <- ggplot(emb, aes(PC1, PC2, col=Microenv))+
  geom_point(size=3)+
  theme_light()+
  scale_color_distiller(palette='Blues',direction = 1)


pcam_h <- ggplot(emb, aes(PC1, PC2, col=BC16_Tumor))+
  geom_point(size=3)+
  theme_light()+
  scale_color_distiller(palette='Reds',direction = 1)

pcat_h <- ggplot(emb, aes(PC1, PC2, col=BC16_Microenv))+
  geom_point(size=3)+
  theme_light()+
  scale_color_distiller(palette='Blues',direction = 1)


pcat/pcam
pcat_h / pcam_h
pcat + pcam + pcat_h + pcam_h



#fix metastasis
clin$MetastasisAtDiagnosis <- "Metastatic"
clin[grepl('Non', clin$`Disease at diagnosis`),"MetastasisAtDiagnosis"] <- "Non-metastatic"





#module score 
markerfile <- '/Users/ferrenaa/Dropbox/data/bangdata/scrnaseq-TKO-DKOAA-DKO/results/Sample-04_DJ582M11/celltypeid/markers-celltypes_lr_2021.07.22.csv'
m <- read.csv(markerfile)

modlist <- list()
plotlist <- list()
celltypes <- unique(m$cluster)
for(ct in celltypes){
  message('Module score for cluster: ', ct)
  
  #get markers for the celltype
  mt <- m[m$cluster==ct,]
  
  #select genes
  mt <- mt[mt$p_val_adj < 0.05,]
  mt <- mt[mt$avg_log2FC > 0.25,]
  mt <- mt[mt$pct.1 > 0.5,]
  mt <- mt[mt$pct.2 < 0.5,]
  
  
  #select orthologs
  genes <- mt$gene
  hugo <- convertMouseGeneList(genes)
  
  #remove duplicates...
  hugo <- hugo[!duplicated(hugo[,1]),]
  hugo <- hugo[!duplicated(hugo[,2]),]
  
  mod <- FerrenaBulkRNAseq::modulescore(gem, hugo$HGNC.symbol)
  
  #do cor
  ctdot <- gsub('-', '.', ct)
  
  ctfinal <- paste0('modulescore_', ctdot)
  
  
  ctdf <- data.frame(cib = clin[,ctdot], mod = mod)
  
  cres <- cor.test(ctdf[,1], ctdf[,2])
  nicep <- ifelse(test = cres$p.value < 0.01, formatC(cres$p.value, format='e', digits=2), 
                  no = round(cres$p.value, 2))
  
  exactlabels <- c(paste0('CIBERSORTx_', ctdot), 
                   ctfinal)
  
  g <- ggplot(ctdf, aes(cib, mod))+
    geom_point(size=3)+
    labs(title = paste0(exactlabels[1], ' vs ', exactlabels[2]),
         subtitle = paste0('Pearson cor r = ', round(cres$estimate, 2),
                           '\nP (cor.test) = ', nicep))+
    ylab(exactlabels[1])+ xlab(exactlabels[2])+
    theme_linedraw() 
  
  plotlist[[ct]] <- g
  modlist[[ct]] <- mod
  
}

patchwork::wrap_plots(plotlist)



# cibsigfile <- '~/Downloads/CIBERSORTx_Job1_cibersortmat_inferred_phenoclasses.CIBERSORTx_Job1_cibersortmat_inferred_refsample.bm.K999.txt'
# c <- read.table(cibsigfile, header = T)
# 








#test metastasis
time <- "Event Free Survival Time in Days"

#set status
if(time == "Event Free Survival Time in Days"){
  status = clin$`First Event`
} else if(time == "Overall Survival Time in Days" ){
  status = clin$`Vital Status`
}


df <- data.frame(samps = sdf$samps,
                 cases = sdf$cases,
                 MetastasisAtDiagnosis = clin$MetastasisAtDiagnosis,
                 surv = clin[,time],
                 status = status)

df$surv <- df$surv / 30

df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)

fit <- survfit(Surv(surv, status_code) ~ MetastasisAtDiagnosis, data = df)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           linetype = "strata",
           surv.median.line = "hv", 
           ggtheme = theme_light(), 
           palette = RColorBrewer::brewer.pal('Set1',n=3)[1:2])+ xlab('Time (Months)')




### event free surv ###
time <- "Event Free Survival Time in Days"

#set status
if(time == "Event Free Survival Time in Days"){
  status = clin$`First Event`
} else if(time == "Overall Survival Time in Days" ){
  status = clin$`Vital Status`
}

# loop thru each cibersort result #

cibersortscores <- names(clin)[31:54]
cibersortscores <- names(clin)[!names(clin) %in% c(origclin,'MetastasisAtDiagnosis')]






plotlist <- list()
reslist <- list()
for(scorename in cibersortscores){
  
  message(scorename)
  
  scorevec <- clin[,scorename]
  
  highname <- paste0('High\nN=', table(scorevec >= median(scorevec) )[2] )
  lowname <- paste0('Low\nN=', table(scorevec >= median(scorevec) )[1] )
  
  dichscore <- ifelse(scorevec >= median(scorevec), yes = highname, no = lowname)
  dichscore <- factor(dichscore, levels = c(lowname, highname))
  
  
  #first dichotomize the score, then run stratified comparison with log-rank test
  df <- data.frame(samps = sdf$samps,
                   cases = sdf$cases,
                   Cibersortscore = scorevec,
                   DichotomizedScore = dichscore,
                   MetastasisAtDiagnosis = clin$MetastasisAtDiagnosis,
                   surv = clin[,time],
                   status = status)
  
  
  #df$DichotomizedScore <- factor(df$DichotomizedScore, levels = rev(unique(df$DichotomizedScore)))
  
  df$surv <- df$surv / 30
  
  df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)
  
  fit <- survfit(Surv(surv, status_code) ~ DichotomizedScore, data = df)
  
  g1 <- ggsurvplot(fit,
                   pval = TRUE, conf.int = TRUE,
                   linetype = "strata", 
                   subtitle=time,
                   title=scorename,
                   surv.median.line = "hv", 
                   ggtheme = theme_light(), 
                   palette = rev( RColorBrewer::brewer.pal('Set1',n=3)[1:2]) )+ 
    xlab('Time (Months)')
  
  
  
  ### cox modelling ###
  #dich model, will match plot
  dich <- coxph(Surv(surv, status_code) ~ DichotomizedScore, data = df)
  
  #continuous univariate cox with cibersort score
  univar <- coxph(Surv(surv, status_code) ~ Cibersortscore, data = df)
  
  #test if there is a difference in scores between mets or not
  metscores <- df[df$MetastasisAtDiagnosis=='Metastatic',"Cibersortscore"]
  nonscores <- df[df$MetastasisAtDiagnosis!='Metastatic',"Cibersortscore"]
  scoremettest <- t.test(metscores, nonscores)
  
  #multivar cox test
  multivar <- coxph(Surv(surv, status_code) ~ MetastasisAtDiagnosis + Cibersortscore, data = df)
  
  
  
  #summarydf
  resdf <- data.frame(
    CIBERSORTxScore = scorename,
    dichHR = exp(dich$coefficients), dichWald = summary(dich)$coefficients[,5],
    univarHR = exp(univar$coefficients), univarWald = summary(univar)$coefficients[,5],
    scoreMetTest = scoremettest$p.value,
    multivarHR = exp( multivar$coefficients[2] ), multivarwald = summary(multivar)$coefficients[2,5] 
  )
  
  plotlist[[scorename]] <- g1
  reslist[[scorename]] <- resdf
  
  
}


res <- dplyr::bind_rows(reslist)

rownames(res) <- NULL

res$significant <- ''
res[res$dichWald<0.05,"significant"] <- '*'
res[res$dichWald<0.001,"significant"] <- '***'

sigres <- res[res$dichWald<0.05,]
sigres

write.csv(res, 'data/survivalresults/eventfreesurvival.csv', row.names = F, quote = F)

pdf('data/survivalresults/eventfreesurvival.pdf')
plotlist
dev.off()







### event free surv ###
time <- "Overall Survival Time in Days"

#set status
if(time == "Event Free Survival Time in Days"){
  status = clin$`First Event`
} else if (time == "Overall Survival Time in Days" ){
  status = clin$`Vital Status`
}

# loop thru each cibersort result #

cibersortscores <- names(clin)[31:54]
cibersortscores <- names(clin)[!names(clin) %in% c(origclin,'MetastasisAtDiagnosis')]






plotlist <- list()
reslist <- list()
for(scorename in cibersortscores){
  
  message(scorename)
  
  scorevec <- clin[,scorename]
  
  highname <- paste0('High\nN=', table(scorevec >= median(scorevec) )[2] )
  lowname <- paste0('Low\nN=', table(scorevec >= median(scorevec) )[1] )
  
  dichscore <- ifelse(scorevec >= median(scorevec), yes = highname, no = lowname)
  dichscore <- factor(dichscore, levels = c(lowname, highname))
  
  
  #first dichotomize the score, then run stratified comparison with log-rank test
  df <- data.frame(samps = sdf$samps,
                   cases = sdf$cases,
                   Cibersortscore = scorevec,
                   DichotomizedScore = dichscore,
                   MetastasisAtDiagnosis = clin$MetastasisAtDiagnosis,
                   surv = clin[,time],
                   status = status)
  
  
  #df$DichotomizedScore <- factor(df$DichotomizedScore, levels = rev(unique(df$DichotomizedScore)))
  
  df$surv <- df$surv / 30
  
  df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)
  
  fit <- survfit(Surv(surv, status_code) ~ DichotomizedScore, data = df)
  
  g1 <- ggsurvplot(fit,
                   pval = TRUE, conf.int = TRUE,
                   linetype = "strata", 
                   subtitle=time,
                   title=scorename,
                   surv.median.line = "hv", 
                   ggtheme = theme_light(), 
                   palette = rev( RColorBrewer::brewer.pal('Set1',n=3)[1:2]) )+ 
    xlab('Time (Months)')
  
  
  
  ### cox modelling ###
  #dich model, will match plot
  dich <- coxph(Surv(surv, status_code) ~ DichotomizedScore, data = df)
  
  #continuous univariate cox with cibersort score
  univar <- coxph(Surv(surv, status_code) ~ Cibersortscore, data = df)
  
  #test if there is a difference in scores between mets or not
  metscores <- df[df$MetastasisAtDiagnosis=='Metastatic',"Cibersortscore"]
  nonscores <- df[df$MetastasisAtDiagnosis!='Metastatic',"Cibersortscore"]
  scoremettest <- t.test(metscores, nonscores)
  
  #multivar cox test
  multivar <- coxph(Surv(surv, status_code) ~ MetastasisAtDiagnosis + Cibersortscore, data = df)
  
  
  
  #summarydf
  resdf <- data.frame(
    CIBERSORTxScore = scorename,
    dichHR = exp(dich$coefficients), dichWald = summary(dich)$coefficients[,5],
    univarHR = exp(univar$coefficients), univarWald = summary(univar)$coefficients[,5],
    scoreMetTest = scoremettest$p.value,
    multivarHR = exp( multivar$coefficients[2] ), multivarwald = summary(multivar)$coefficients[2,5] 
  )
  
  plotlist[[scorename]] <- g1
  reslist[[scorename]] <- resdf
  
  
}


res <- dplyr::bind_rows(reslist)

rownames(res) <- NULL

res$significant <- ''
res[res$dichWald<0.05,"significant"] <- '*'
res[res$dichWald<0.001,"significant"] <- '***'

sigres <- res[res$dichWald<0.05,]
sigres

write.csv(res, 'data/survivalresults/overallsurvival.csv', row.names = F, quote = F)

pdf('data/survivalresults/overallsurvival.pdf')
plotlist
dev.off()



### downstream analysis ###




plotlist[['X5']] 
plotlist[['Monocyte.Macrophage']]
plotlist[['MSC.Fibroblast']]

plotlist[['X1']]
plotlist[['X7']]
plotlist[['X4']]

plotlist[['BC16_Macrophage']]
plotlist[['BC16_Fibroblast']]




qplot(clin$MetastasisAtDiagnosis, clin$Monocyte.Macrophage, geom='boxplot')
qplot(clin$MetastasisAtDiagnosis, clin$X5, geom='boxplot')
qplot(clin$MetastasisAtDiagnosis, clin$MSC.Fibroblast, geom='boxplot')
qplot(clin$MetastasisAtDiagnosis, clin$X4, geom='boxplot')




scorename='X4'

message(scorename)

scorevec <- clin[,scorename]

highname <- paste0('High\nN=', table(scorevec >= median(scorevec) )[2] )
lowname <- paste0('Low\nN=', table(scorevec >= median(scorevec) )[1] )

dichscore <- ifelse(scorevec >= median(scorevec), yes = highname, no = lowname)
dichscore <- factor(dichscore, levels = c(lowname, highname))


#first dichotomize the score, then run stratified comparison with log-rank test
df <- data.frame(samps = sdf$samps,
                 cases = sdf$cases,
                 Cibersortscore = scorevec,
                 DichotomizedScore = dichscore,
                 MetastasisAtDiagnosis = clin$MetastasisAtDiagnosis,
                 surv = clin[,time],
                 status = clin$`Vital Status`)

df = df[df$MetastasisAtDiagnosis=='Metastatic',]

#df$DichotomizedScore <- factor(df$DichotomizedScore, levels = rev(unique(df$DichotomizedScore)))

df$surv <- df$surv / 30

df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)

fit <- survfit(Surv(surv, status_code) ~ DichotomizedScore, data = df)

g1 <- ggsurvplot(fit,
                 pval = TRUE, conf.int = TRUE,
                 linetype = "strata", 
                 subtitle=time,
                 title=scorename,
                 surv.median.line = "hv", 
                 ggtheme = theme_light(), 
                 palette = rev( RColorBrewer::brewer.pal('Set1',n=3)[1:2]) )+ 
  xlab('Time (Months)')



#first dichotomize the score, then run stratified comparison with log-rank test
df <- data.frame(samps = sdf$samps,
                 cases = sdf$cases,
                 Cibersortscore = scorevec,
                 DichotomizedScore = dichscore,
                 MetastasisAtDiagnosis = clin$MetastasisAtDiagnosis,
                 surv = clin[,time],
                 status = clin$`Vital Status`)

df = df[df$MetastasisAtDiagnosis!='Metastatic',]

#df$DichotomizedScore <- factor(df$DichotomizedScore, levels = rev(unique(df$DichotomizedScore)))

df$surv <- df$surv / 30

df$status_code <- ifelse(df$status == 'Dead', yes=2, no = 1)

fit <- survfit(Surv(surv, status_code) ~ DichotomizedScore, data = df)

g2 <- ggsurvplot(fit,
                 pval = TRUE, conf.int = TRUE,
                 linetype = "strata", 
                 subtitle=time,
                 title=scorename,
                 surv.median.line = "hv", 
                 ggtheme = theme_light(), 
                 palette = rev( RColorBrewer::brewer.pal('Set1',n=3)[1:2]) )+ 
  xlab('Time (Months)')



#try to look at histoligic response and cluster 1

clinx <- clin[,c(26,27,41)]
head(clinx)


clinx <- clinx[!is.na(clinx$`Histologic response`) | !is.na(clinx$`Percent necrosis at Definitive Surgery`),]


#fix perc necrosis
clinx[6, 2] <- 50
clinx[,2] <- as.numeric(clinx$`Percent necrosis at Definitive Surgery`)


#fix histologic response
clinx[grepl('Stage', clinx$`Histologic response`),"Histologic response"] <- NA



clinx$responder <- NA
nonresponder <- ifelse( na.omit(clinx$`Percent necrosis at Definitive Surgery`) >= 91, 
                        yes = "responder", no = 'non-responder')

nonresponder <- c(nonresponder,
                  ifelse(na.omit(clinx$`Histologic response`) == '91-100',
                         yes = 'responder', no = 'non-responder')
                  )


clinx$responder <- nonresponder





