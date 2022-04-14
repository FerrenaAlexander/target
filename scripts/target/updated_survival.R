library(tidyverse)
library(preprocessCore)
library(survival)
library(survminer)
library(DESeq2)
library(patchwork)

library(biomaRt) #mouse-human homologs

library(MCPcounter) #for calculating microenv scores in human

set.seed(2020)




#set up biomart function for mouse to human
human = NULL; mouse = NULL

while(is(human) != 'Mart'){ human = try( useMart("ensembl", dataset = "hsapiens_gene_ensembl")) }

while(is(mouse) != 'Mart'){ mouse = try( useMart("ensembl", dataset = "mmusculus_gene_ensembl")) }




convertMouseGeneList <- function(genelist){
  require("biomaRt")
  
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), 
                   filters = "mgi_symbol", 
                   values = genelist , 
                   mart = mouse, 
                   attributesL = c("hgnc_symbol"), 
                   martL = human, 
                   uniqueRows=T)
  
  genesV2 
}



# set up module score function

modulescore <- function(gem, genelist, numbins=NULL, numcontrolgenesperbin=NULL){
  if(is.null(numbins)){numbins=24}
  if(is.null(numcontrolgenesperbin)){numcontrolgenesperbin=100}
  
  testgenes <- genelist
  controlgenes <- rownames(gem)
  controlgenes <- controlgenes[!(controlgenes %in% testgenes)]
  
  #testmeans <- Matrix::rowMeans(gem[rownames(gem) %in% genes,])
  
  ### bin the genes by avg expression ###
  
  #get average gene expression
  avgs <- Matrix::rowMeans(gem)
  avgs <- avgs[order(avgs)]
  
  #cut; get the gene bins
  bins <- cut_number(avgs, n=numbins)
  
  bindf <- data.frame(genenames = names(avgs), 
                      bins = bins)
  
  #select control genes from same expression bins
  controlgenes <- c()
  for(genedex in 1:length(testgenes)){
    
    #get gene and bin
    gene <- testgenes[genedex]
    genebin <- bindf[bindf$genenames == gene,2]
    
    #select all genes from that bin, to sample from
    tmpcontrols <- bindf[bindf$bins == genebin,]
    
    #exclude the actual test genes
    tmpcontrols <- tmpcontrols[!(tmpcontrols$genenames %in% testgenes),]
    
    #if num controls exceeds bin size, take all genes in bin...
    numtotake <- ifelse(nrow(tmpcontrols) < numcontrolgenesperbin,
                        yes = nrow(tmpcontrols),
                        no = numcontrolgenesperbin)
    
    #sample the genes
    tmpcontrols <- sample(tmpcontrols$genenames, size = numtotake, replace = F)
    
    controlgenes <- unique(c(controlgenes, 
                             tmpcontrols
    ))
  }
  
  
  
  #get control gene mean expression for each sample
  samplemeans_controls <- Matrix::colMeans( gem[rownames(gem) %in% controlgenes,] )
  
  #get test genes mean expression for each samoke
  samplemeans_test <- Matrix::colMeans( gem[rownames(gem) %in% testgenes,] )
  
  #subtract them to get module score
  modulescore <- samplemeans_test - samplemeans_controls
  
  return(modulescore)
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

rm(l, sampmat)



#fix clin

# binarize First Event status to censored / event (ie no / yes or 0/1)
firstevent <- clin$`First Event`
firstevent[is.na(firstevent)] <- 'NA'
firstevent[ firstevent == 'Censored' | firstevent == 'None' ] <- 'Censored'
firstevent[ firstevent!= 'NA' & firstevent != 'Censored' ] <- 'Event'
firstevent[firstevent == 'NA'] <- NA

clin$event_recoded <- factor(firstevent, levels = c('Event', 'Censored') )
rm(firstevent)

# recode time to first event
firsteventtime <- clin$`Event Free Survival Time in Days`
firsteventtime <- firsteventtime / 365.25

clin$event_recoded_time <- firsteventtime 
rm(firsteventtime)

# recode vital status to death 1 = yes, 0 = no
vitalstat <- clin$`Vital Status`
vitalstat[is.na(vitalstat)] <- 'NA'
vitalstat[ vitalstat == 'Alive'] <- 'Censored'
vitalstat[ vitalstat!= 'NA' & vitalstat != 'Censored' ] <- 'Dead'
vitalstat[vitalstat == 'NA'] <- NA

clin$vitalstatus_recoded <- factor(vitalstat, levels = c('Dead', 'Censored') )
rm(vitalstat)

# recode survival time
ostime <- clin$`Overall Survival Time in Days`
ostime <- ostime / 365.25

clin$vitalstatus_recoded_time <- ostime 
rm(ostime)


# recode metastasis, remove parenthesis etc
clin$MetastasisAtDiagnosis <- "Metastatic"
clin[grepl('Non', clin$`Disease at diagnosis`),"MetastasisAtDiagnosis"] <- "Non-metastatic"

clin$MetastasisAtDiagnosis <- factor(clin$MetastasisAtDiagnosis, levels = c('Non-metastatic', 'Metastatic'))

#recode age to years
clin$age_years <- clin$`Age at Diagnosis in Days` / 365.25



#set up cleaner clinical variable dataframe

clindat <- data.frame(Sample = clin$`TARGET USI`,
                      Gender = clin$Gender,
                      Age_at_diagnosis = clin$age_years,
                      
                      Mestastasis_at_diagnosis = clin$MetastasisAtDiagnosis,
                      
                      First_event = clin$event_recoded,
                      First_event_time = clin$event_recoded_time,
                      
                      Vital_status = clin$vitalstatus_recoded,
                      Vital_status_time = clin$vitalstatus_recoded_time
                      )








rm(clin)


#check NAs
clindat[!complete.cases(clindat),]

#one sample can be fixed: TARGET-40-0A4I9K
# marked as dead, but event is missing
# even tho event is missing, first event time is there
# vital status time and first event time are both >= 5
# we will set first event as true for this sample.
clindat[clindat$Sample == 'TARGET-40-0A4I9K',"First_event"] <- "Event"

#other samples are missing any survival, just remove them.

clindat <- clindat[complete.cases(clindat),]



#cap at 5 years
# set time above 5 to just 5
# also we need to censor after 5 years before doing that

clindat[clindat$First_event_time >5,'First_event'] <- 'Censored'
clindat[clindat$First_event_time >5,'First_event_time'] <- 5


clindat[clindat$Vital_status_time >5,'Vital_status'] <- 'Censored'
clindat[clindat$Vital_status_time >5,'Vital_status_time'] <- 5



#remove samples missing survival

sdf <- sdf[match(clindat$Sample, sdf$cases),]
gem <- gem[,match(sdf$samps, colnames(gem))]

# also rename gem colnames... from "sample" to "case"
colnames(gem) <- clindat$Sample
rm(sdf)


#check lib size
libsize <- data.frame(samp = colnames(gem),
                      libsize = colSums(gem))


ggplot(libsize, aes(x = libsize))+
  geom_histogram(col = 'black', fill='steelblue')+
  scale_x_log10()



#remove extreme low TPM samples
badsamps <- c('TARGET-40-PALFYN',
              'TARGET-40-PASKZZ')

clindat <- clindat[!( clindat$Sample %in% badsamps),]
gem <- gem[,colnames(gem) %in% clindat$Sample]

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






#### calculate scores ####

# do this before PCA, we want to see how they loook in PCA



#calculate the module scores

#read in DE genes
res <- read.csv('~/Dropbox/data/bangdata/2021october-TKOvsDKO/results/comparative-de/TKO-vs-DKO/2.deresults/DEresults-normalized_count_matrix.csv')

#get only sig res
res <- res[abs(res$log2FoldChange) > 1 & res$padj < 0.05,]


#get the human IDs...

#remove dup genes (if any)
dups <- names( table(res$mgi_symbol)[table(res$mgi_symbol)>1] )
res <- res[!(res$mgi_symbol %in% dups),]


#get orthologs
genes <- convertMouseGeneList(res$mgi_symbol)

#check dups...
table(duplicated(genes$MGI.symbol))
table(duplicated(genes$HGNC.symbol))

dup_mgi <- names( table(genes$MGI.symbol)[table(genes$MGI.symbol)>1] )
genes[genes$MGI.symbol %in% dup_mgi,]

dup_hgnc <- names( table(genes$HGNC.symbol)[table(genes$HGNC.symbol)>1] )
genes[genes$HGNC.symbol %in% dup_hgnc,]

#can't figure out how to resolve dups, just remove for now, not ideal...
genes_nodups <- genes[!(genes$MGI.symbol %in% dup_mgi),]
genes_nodups <- genes_nodups[!(genes_nodups$HGNC.symbol %in% dup_hgnc),]



#try to match up...
res <- res[res$mgi_symbol %in% genes_nodups$MGI.symbol,]
genes_nodups <- genes_nodups[match(res$mgi_symbol, genes_nodups$MGI.symbol),]
res <- cbind(genes_nodups$HGNC.symbol, res)
colnames(res)[1] <- 'HGNC_ortholog'

#split into up/down lists
up <- res[res$log2FoldChange>0,"HGNC_ortholog"]
down <- res[res$log2FoldChange<0,"HGNC_ortholog"]

mm <- data.frame(TKO_overexpressed = modulescore(gem, up),
                 TKO_underexpressed = modulescore(gem, down))






#MCP counter scores

# links to a live github: 2022.01.05
# mcp <- MCPcounter.estimate(gem, featuresType = 'HUGO_symbols')
# saveRDS(mcp, 'data/mcpresults_2022.rds')
mcp <- readRDS("data/mcpresults_2022.rds")

#check mcp gene markers...
# links to a live github: 2022.01.05
#genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)



#add to mm
rownames(mcp) <- gsub(' ', '_', rownames(mcp))
mm <- cbind(mm,t(mcp))


#plot it...

mcpplot<- t(mcp)
mcpplot <- reshape2::melt(mcpplot)

colnames(mcpplot) <- c('Sample', 'Celltype', 'Estimate')

g1 <- ggplot(mcpplot, aes(Celltype, Estimate))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

g2 <- ggplot(mcpplot, aes(Celltype, Estimate))+
  geom_boxplot()+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

g1 / g2


#exclude fibroblasts...
mcpplot<- t(mcp)
mcpplot <- reshape2::melt(mcpplot[,-10])

colnames(mcpplot) <- c('Sample', 'Celltype', 'Estimate')

g1 <- ggplot(mcpplot, aes(Celltype, Estimate))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

g2 <- ggplot(mcpplot, aes(Celltype, Estimate))+
  geom_boxplot()+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

g1 / g2



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

emb <- emb[,1:2]

#get genes
pcagenes <- pca$rotation
head(sort(pcagenes[,1],decreasing = T))
head(sort(pcagenes[,2],decreasing = T))

libsize <- libsize[libsize$samp %in% clindat$Sample,]
emb$libsize <- libsize$libsize

#add other clinical vars
emb <- cbind(emb, clindat)

#add module scores
emb <- cbind(emb, scale(mm))

#plot PCs
# with lib size
ggplot(emb, aes(PC1, PC2, col=libsize))+
  geom_point()+
  theme_light()+
  scale_color_distiller(palette = 'Purples')

# by sex
ggplot(emb, aes(PC1, PC2, col=Gender))+
  geom_point()+
  theme_light()+
  scale_color_brewer(palette = 'Set1', direction = 1)

# by age
ggplot(emb, aes(PC1, PC2, col=Age_at_diagnosis))+
  geom_point()+
  theme_light()+
  scale_color_distiller(palette = 'Greens', direction = 1)


# by stage; ie met or not
ggplot(emb, aes(PC1, PC2, col=Mestastasis_at_diagnosis))+
  geom_point()+
  theme_light()+
  scale_color_brewer(palette = 'Set2')


# by event status
ggplot(emb, aes(PC1, PC2, col=First_event))+
  geom_point()+
  theme_light()+
  scale_color_brewer(palette = 'Set2') +
  
  ggplot(emb, aes(PC1, PC2, col=First_event_time))+
  geom_point()+
  theme_light()+
  scale_color_distiller(palette = 'Purples')


# by vital status
ggplot(emb, aes(PC1, PC2, col=Vital_status))+
  geom_point()+
  theme_light()+
  scale_color_brewer(palette = 'Set2') +
  
  ggplot(emb, aes(PC1, PC2, col=Vital_status_time))+
  geom_point()+
  theme_light()+
  scale_color_distiller(palette = 'Purples')


#by module scores
ggplot(emb, aes(PC1, PC2, col=TKO_overexpressed))+
  geom_point()+
  theme_light()+
  scale_color_distiller(palette = 'Reds') +
  
  ggplot(emb, aes(PC1, PC2, col=TKO_underexpressed))+
  geom_point()+
  theme_light()+
  scale_color_distiller(palette = 'Blues')


# cell types...
ct <- rownames(mcp)
ctplots <- list()
for(t in ct){
  df <- data.frame(PC1 = emb$PC1,
                   PC2 = emb$PC2,
                   celltype = emb[,t])
  ctplots[[t]] <- ggplot(df, aes(PC1, PC2, col=celltype))+
    geom_point()+
    theme_light()+
    scale_color_distiller(palette = 'RdBu', name = t)
  rm(df)
}
patchwork::wrap_plots(ctplots)




rm(emb, gem_hvgs, gem_scale, gemt, libsize_post, pcagenes, pca, numhvgs, rowvars)


#lets use vital status
status <- 'Vital_status'
time <- 'Vital_status_time'



#set outdir
outdir <- 'survivalanalysis/vital_status'
dir.create(outdir)


### survival analysis for clinical variables


suboutdir <- paste0(outdir, '/clinicalvars')
dir.create(suboutdir)

# start with continuous vars
varnames <- names(clindat)
varnames <- varnames[c(3)]
#age is the only one? lol


plotlist <- list()
modellist <- list()
for(scorename in varnames){
  
  
  message('\n', scorename)
  
  scorevec <- clindat[,scorename]
  
  highname <- paste0('High\nN=', table(na.omit(scorevec) >= median(na.omit(scorevec)) )[2] )
  lowname <- paste0('Low\nN=', table(na.omit(scorevec) >= median(na.omit(scorevec)) )[1] )
  
  dichscore <- ifelse(scorevec >= median(na.omit(scorevec)), yes = highname, no = lowname)
  dichscore <- factor(dichscore, levels = c(lowname, highname))
  
  
  
  #add terciles
  terciles <- factor(ntile(scorevec, 3))
  terciles <- plyr::mapvalues(x = terciles, from = levels(terciles),
                              to = c( paste0('Low\nN=', table(terciles)[1]),
                                      paste0('Med\nN=', table(terciles)[2]),
                                      paste0('High\nN=', table(terciles)[2]) )
  )
  
  
  
  df <- data.frame(samps = clindat$Sample,
                   continuous = scorevec,
                   dichotomized = dichscore,
                   terciles = terciles,
                   surv = clindat[,time],
                   status = clindat[,status])
  
  
  #numeric cding for vital status
  df$status_code <- ifelse(df$status == 'Censored', yes=0, no = 1)
  
  
  fit <- survfit(Surv(surv, status_code) ~ terciles, data = df)
  
  gterc <- ggsurvplot(fit,
                      conf.int = F,
                      pval = TRUE, pval.coord = c(3, 0.9),
                      pval.method = T, pval.method.coord = c(3, 1),
                      linetype = "strata",  
                      test.for.trend = T,
                      title=scorename,
                      ggtheme = theme_classic2(), 
                      palette =  c( "#8DA0CB", "#66C2A5", "#FC8D62" ) )+ 
    xlab('Time (Years)')
  
  fit <- survfit(Surv(surv, status_code) ~ dichotomized, data = df)
  
  gdich <- ggsurvplot(fit,
                      conf.int = F,
                      pval = TRUE, pval.coord = c(3, 0.9),
                      pval.method = T, pval.method.coord = c(3, 1),
                      linetype = "strata",  
                      test.for.trend = F,
                      title=scorename,
                      ggtheme = theme_classic2(), 
                      palette = c("#377EB8" , "#E41A1C") )+ 
    xlab('Time (Years)')
  
  
  gterc <- gterc$plot
  gdich <- gdich$plot
  
  ### cox modelling ###
  #dich model, will match plot
  dich <- coxph(Surv(surv, status_code) ~ dichotomized, data = df)
  
  #terc
  terc <- coxph(Surv(surv, status_code) ~ terciles, data = df)
  
  #continuous univariate cox with cibersort score
  continuous <- coxph(Surv(surv, status_code) ~ continuous, data = df)
  
  #multivar cox test --> not for clinical
  # multivar <- coxph(Surv(surv, status_code) ~ stage + age, data = df)
  
  #save all outputs...
  resultdir <- paste0(suboutdir, '/', scorename)
  dir.create(resultdir)
  
  #save models
  models <- list(dich = dich,
                 terc = terc,
                 continuous = continuous
                 #multivar = multivar
  )
  
  saveRDS( file = paste0(resultdir, '/models.rds'), models  )
  
  
  
  #check model assumptions...
  modeldiagnosticsfile <- paste0(resultdir, '/modeldiagnostics.pdf')
  pdf(modeldiagnosticsfile, height = 10, width = 10)
  for(model in models){
    
    print( ggcoxzph( cox.zph(model) ) )
    
  }
  dev.off()
  
  
  
  
  #save plots
  plots <- list(gdich, gterc)
  saveRDS( file = paste0(resultdir, '/plots.rds'), plots   )
  
  ggsave( paste0(resultdir, '/gdich.jpg'), gdich , height = 7, width = 7)
  ggsave( paste0(resultdir, '/gterc.jpg'), gterc , height = 7, width = 7)
  
  #add all to lists
  plotlist[[scorename]] <- plots
  modellist[[scorename]] <- models
  
  
}





#categorical vars
varnames <- names(clindat)
varnames <- varnames[varnames != 'Sample']
varnames <- varnames[c(1,3)]


for(scorename in varnames){
  
  message('\n', scorename)
  
  scorevec <- clindat[,scorename]
  
  # change name to have N
  scorevec <- factor(scorevec)
  scorevec <- plyr::mapvalues(scorevec, from = levels(scorevec),
                              to = paste0(levels(scorevec), '\nN=', table(scorevec)[levels(scorevec)] ) )
  
  
  
  
  df <- data.frame(samps = clindat$Sample,
                   var = scorevec,
                   surv = clindat[,time],
                   status = clindat[,status])

  # code status as numeric
  df$status_code <- ifelse(df$status == 'Censored', yes=0, no = 1)
  
  fit <- survfit(Surv(surv, status_code) ~ var, data = df)
  
  gcat <- ggsurvplot(fit,
                     conf.int = F,
                     pval = TRUE, pval.coord = c(3, 0.9),
                     pval.method = T, pval.method.coord = c(3, 1),
                     linetype = "strata",  
                     test.for.trend = F,
                     title=scorename,
                     #surv.median.line = "hv", 
                     ggtheme = theme_classic2())+ 
    xlab('Time (Years)')
  
  #test for trend, b ut only if >2 vars...
  if( length(levels(df$var)) > 2 )  {
    
    gtft <- ggsurvplot(fit,
                       conf.int = F,
                       pval = TRUE, pval.coord = c(3, 0.9),
                       pval.method = T, pval.method.coord = c(3, 1),
                       linetype = "strata",  
                       test.for.trend = T,
                       title=scorename,
                       #surv.median.line = "hv", 
                       ggtheme = theme_classic2())+ 
      xlab('Time (Years)')
    
  }
  
  gcat <- gcat$plot
  
  if( length(levels(df$var)) > 2 )  {
    gtft <- gtft$plot
    plots <- list(gcat, gtft)
    
  } else{
    plots <- gcat
  }
  
  ### cox modelling ###
  #dich model, will match plot
  model <- coxph(Surv(surv, status_code) ~ var, data = df)
  
  
  #just save the likelihood p value and assess the saved model
  resdf <- data.frame(
    Module = scorename,
    modelP = summary(model)$logtest[3],
    #multivarP = summary(multivar)$logtest[3],
    row.names = NULL
    
  )
  
  
  #save all outputs...
  resultdir <- paste0(suboutdir, '/', scorename)
  dir.create(resultdir)
  
  #save models -> just one for categoricals
  models <- list(model = model
                 
  )
  
  
  #check model assumptions...
  modeldiagnosticsfile <- paste0(resultdir, '/modeldiagnostics.pdf')
  pdf(modeldiagnosticsfile, height = 10, width = 10)
  for(model in models){
    
    print( ggcoxzph( cox.zph(model) ) )
    
  }
  dev.off()
  
  
  
  saveRDS( file = paste0(resultdir, '/models.rds'), models  )
  
  #save plots
  saveRDS( file = paste0(resultdir, '/plots.rds'), plots   )
  
  ggsave( paste0(resultdir, '/gdich.jpg'), gdich , height = 7, width = 7)
  ggsave( paste0(resultdir, '/gterc.jpg'), gterc , height = 7, width = 7)
  
  #add all to lists
  plotlist[[scorename]] <- plots
  modellist[[scorename]] <- model
  
}



#summary

summaryfile <- paste0(suboutdir, '/summary-plots_models.pdf')
pdf(file = summaryfile, height = 10,width = 9)

for(score in names(plotlist) ){
  plots <- plotlist[[score]]
  
  #get models liek this because continous have 3 cats have 1
  models <- modellist[score]
  
  #continus have 3 models, cats have 1
  if( length(models[[1]]) == 3 ){
    
    models <- modellist[[score]]
    summarytable <- bind_rows( lapply(models, function(x){as.data.frame( summary(x)$coefficients ) }) )
    summarytable <- bind_cols(summarytable,
                              bind_rows( lapply(models, function(x){as.data.frame( summary(x)$conf.int ) }) )[,3:4] )
    
  } else{
    summarytable <- bind_rows( lapply(models, function(x){as.data.frame( summary(x)$coefficients ) }) )
    summarytable <- bind_cols(summarytable,
                              bind_rows( lapply(models, function(x){as.data.frame( summary(x)$conf.int ) }) )[,3:4] )
    
  }
  
  #put conf int next to coef
  summarytable <- summarytable[,c(1,6,7,2,3,4,5)]
  
  
  #add significance marks
  summarytable$significance <- ''
  summarytable[summarytable$`Pr(>|z|)` < 0.1, "significance"] <- '.'
  summarytable[summarytable$`Pr(>|z|)` < 0.05, "significance"] <- '*'
  summarytable[summarytable$`Pr(>|z|)` < 0.01, "significance"] <- '**'
  summarytable[summarytable$`Pr(>|z|)` < 0.001, "significance"] <- '***'
  
  #round table
  summarytable[,-ncol(summarytable)] <- round( summarytable[,-ncol(summarytable)],2 )
  
  #osme categorical vars have only one plot
  if(length(plots) ==2 ) {
    updatedplot <- (plots[[1]] + plots[[2]]) / ( gridExtra::tableGrob(summarytable)) + patchwork::plot_annotation(score)
  } else{
    updatedplot <- (plots) / ( gridExtra::tableGrob(summarytable)) + patchwork::plot_annotation(score)
    
  }
  
  print(updatedplot)
  
}
dev.off()


















####################################### module scores survival correlations #############################################


suboutdir <- paste0(outdir, '/modulescores')
dir.create(suboutdir)

#get module score names
scorenames <- names(mm)

plotlist <- list()
modellist <- list()

for(scorename in scorenames){
  
  message(scorename)
  
  scorevec <- mm[,scorename]
  
  
  highname <- paste0('High\nN=', table(na.omit(scorevec) >= median(na.omit(scorevec)) )[2] )
  lowname <- paste0('Low\nN=', table(na.omit(scorevec) >= median(na.omit(scorevec)) )[1] )
  
  dichscore <- ifelse(scorevec >= median(na.omit(scorevec)), yes = highname, no = lowname)
  dichscore <- factor(dichscore, levels = c(lowname, highname))
  
  
  # highname <- paste0('High\nN=', table(na.omit(scorevec) >= 0 )[2] )
  # lowname <- paste0('High\nN=', table(na.omit(scorevec) >= 0 )[2] )
  
  

  
  
  #add terciles
  terciles <- factor(ntile(scorevec, 3))
  terciles <- plyr::mapvalues(x = terciles, from = levels(terciles),
                              to = c( paste0('Low\nN=', table(terciles)[1]),
                                      paste0('Med\nN=', table(terciles)[2]),
                                      paste0('High\nN=', table(terciles)[2]) )
  )
  
  
  

  df <- data.frame(samps = clindat$Sample,
                   continuous = scorevec,
                   dichotomized = dichscore,
                   terciles = terciles,
                   surv = clindat[,time],
                   status = clindat[,status])
  
  
  # code status as numeric
  df$status_code <- ifelse(df$status == 'Censored', yes=0, no = 1)
  
  fit <- survfit(Surv(surv, status_code) ~ terciles, data = df)
  
  gterc <- ggsurvplot(fit,
                      conf.int = F,
                      pval = TRUE, pval.coord = c(3, 0.9),
                      pval.method = T, pval.method.coord = c(3, 1),
                      linetype = "strata",  
                      test.for.trend = T,
                      title=scorename,
                      ggtheme = theme_classic2(), 
                      palette =  c( "#8DA0CB", "#66C2A5", "#FC8D62" ) )+ 
    xlab('Time (Years)')
  
  
  fit <- survfit(Surv(surv, status_code) ~ dichotomized, data = df)
  
  gdich <- ggsurvplot(fit,
                      conf.int = F,
                      pval = TRUE, pval.coord = c(3, 0.9),
                      pval.method = T, pval.method.coord = c(3, 1),
                      linetype = "strata",  
                      test.for.trend = F,
                      title=scorename,
                      ggtheme = theme_classic2(), 
                      palette = c("#377EB8" , "#E41A1C") )+ 
    xlab('Time (Years)')
  
  
  gterc <- gterc$plot
  gdich <- gdich$plot
  
  ### cox modelling ###
  #dich model, will match plot
  dich <- coxph(Surv(surv, status_code) ~ dichotomized, data = df)
  
  #terc
  terc <- coxph(Surv(surv, status_code) ~ terciles, data = df)
  
  #continuous univariate cox with cibersort score
  continuous <- coxph(Surv(surv, status_code) ~ continuous, data = df)
  
  #prep for multivar
  df <- cbind(df, clindat)
  
  #multivar cox test 
  multivar <- coxph(Surv(surv, status_code) ~ continuous + Mestastasis_at_diagnosis, data = df)
  
  
  #save all outputs...
  resultdir <- paste0(suboutdir, '/', scorename)
  dir.create(resultdir)
  
  #save models
  models <- list(dich = dich,
                 terc = terc,
                 continuous = continuous,
                 multivar = multivar
  )
  
  saveRDS( file = paste0(resultdir, '/models.rds'), models  )
  
  
  #check model assumptions...
  modeldiagnosticsfile <- paste0(resultdir, '/modeldiagnostics.pdf')
  pdf(modeldiagnosticsfile, height = 10, width = 10)
  for(model in models){
    
    print( ggcoxzph( cox.zph(model) ) )
    
  }
  dev.off()
  
  
  #save plots
  plots <- list(gdich, gterc)
  saveRDS( file = paste0(resultdir, '/plots.rds'), plots   )
  
  ggsave( paste0(resultdir, '/gdich.jpg'), gdich , height = 7, width = 7)
  ggsave( paste0(resultdir, '/gterc.jpg'), gterc , height = 7, width = 7)
  
  #add all to lists
  plotlist[[scorename]] <- plots
  modellist[[scorename]] <- models
  
}




#summary

summaryfile <- paste0(suboutdir, '/summary-plots_models.pdf')
pdf(file = summaryfile, height = 10,width = 9)

for(score in names(plotlist) ){
  plots <- plotlist[[score]]
  
  #get models
  models <- modellist[[score]]
  
  #get the coefficients, hazard ratios, confidence intervals, and P values
  summarytable <- bind_rows( lapply(models, function(x){as.data.frame( summary(x)$coefficients ) }) )
  summarytable <- bind_cols(summarytable,
                            bind_rows( lapply(models, function(x){as.data.frame( summary(x)$conf.int ) }) )[,3:4] )
  
  
  #rename the continious...
  rownames <- rownames(summarytable)
  rownames[grepl('^continuous', rownames)] <- c('Univariable\nContinuous', 'Multivariable\nContinuous')
  rownames(summarytable) <- rownames
  
  #put conf int next to coef
  summarytable <- summarytable[,c(1,6,7,2,3,4,5)]
  
  
  #add significance marks
  summarytable$significance <- ''
  summarytable[summarytable$`Pr(>|z|)` < 0.1, "significance"] <- '.'
  summarytable[summarytable$`Pr(>|z|)` < 0.05, "significance"] <- '*'
  summarytable[summarytable$`Pr(>|z|)` < 0.01, "significance"] <- '**'
  summarytable[summarytable$`Pr(>|z|)` < 0.001, "significance"] <- '***'
  
  #round table
  summarytable[,-ncol(summarytable)] <- round( summarytable[,-ncol(summarytable)],2 )
  
  # plot it all together
  updatedplot <- (plots[[1]] + plots[[2]]) / ( gridExtra::tableGrob(summarytable)) + patchwork::plot_annotation(score) + 
    plot_layout(heights = c(1,2))
  
  
  print(updatedplot)
  
}
dev.off()















############# stratify by metastasis state ############



suboutdir <- paste0(outdir, '/modulescores_dichotomized_by_metastasis')
dir.create(suboutdir)

#get module score names
scorenames <- names(mm)

outputlist <- list()

for(scorename in scorenames){
  
  message(scorename)
  scorevec <- mm[,scorename]
  
  
  #first, dichotomize by met
  df <- data.frame(samps = clindat$Sample,
                   continuous = scorevec,
                   surv = clindat[,time],
                   status = clindat[,status])
  
  df <- cbind(df, clindat)
  
  #save all outputs
  resultdir <- paste0(suboutdir, '/', scorename)
  dir.create(resultdir)
  
  #run analysis in age strata
  stratumlist <- list()
  for(stratum in levels(df$Mestastasis_at_diagnosis) ) {
    
    sdf <- df[df$Mestastasis_at_diagnosis==stratum, ]
    scorevec <- sdf$continuous
    
    
    
    highname <- paste0('High\nN=', table(na.omit(scorevec) >= median(na.omit(scorevec)) )[2] )
    lowname <- paste0('Low\nN=', table(na.omit(scorevec) >= median(na.omit(scorevec)) )[1] )
    
    dichscore <- ifelse(scorevec >= median(na.omit(scorevec)), yes = highname, no = lowname)
    dichscore <- factor(dichscore, levels = c(lowname, highname))
    
    
    # highname <- paste0('High\nN=', table(na.omit(scorevec) >= 0 )[2] )
    # lowname <- paste0('High\nN=', table(na.omit(scorevec) >= 0 )[2] )
    
    
    
    
    
    #add terciles
    terciles <- factor(ntile(scorevec, 3))
    terciles <- plyr::mapvalues(x = terciles, from = levels(terciles),
                                to = c( paste0('Low\nN=', table(terciles)[1]),
                                        paste0('Med\nN=', table(terciles)[2]),
                                        paste0('High\nN=', table(terciles)[2]) )
    )
    
    
    
    #add to stratum df
    sdf$dichotomized <- dichscore
    sdf$terciles <- terciles
    
    
    
    # code status as numeric
    sdf$status_code <- ifelse(sdf$status == 'Censored', yes=0, no = 1)
    
    fit <- survfit(Surv(surv, status_code) ~ terciles, data = sdf)
    
    gterc <- ggsurvplot(fit,
                        conf.int = F,
                        pval = TRUE, pval.coord = c(3, 0.9),
                        pval.method = T, pval.method.coord = c(3, 1),
                        linetype = "strata",  
                        test.for.trend = T,
                        title=paste0('Stratum: ', stratum),
                        ggtheme = theme_classic2(), 
                        palette =  c( "#8DA0CB", "#66C2A5", "#FC8D62" ) )+ 
      xlab('Time (Years)')
    
    
    fit <- survfit(Surv(surv, status_code) ~ dichotomized, data = sdf)
    
    gdich <- ggsurvplot(fit,
                        conf.int = F,
                        pval = TRUE, pval.coord = c(3, 0.9),
                        pval.method = T, pval.method.coord = c(3, 1),
                        linetype = "strata",  
                        test.for.trend = F,
                        title=paste0('Stratum: ', stratum),
                        ggtheme = theme_classic2(), 
                        palette = c("#377EB8" , "#E41A1C") )+ 
      xlab('Time (Years)')
    
    
    gterc <- gterc$plot
    gdich <- gdich$plot
    
    ### cox modelling ###
    #dich model, will match plot
    dich <- coxph(Surv(surv, status_code) ~ dichotomized, data = sdf)
    
    #terc
    terc <- coxph(Surv(surv, status_code) ~ terciles, data = sdf)
    
    #continuous univariate cox with cibersort score
    continuous <- coxph(Surv(surv, status_code) ~ continuous, data = sdf)
    
    #multivar cox test  --> not needed currently, onyl confounder was met
    # multivar <- coxph(Surv(surv, status_code) ~ continuous + Mestastasis_at_diagnosis, data = df)
    
    #save outputs
    
    #save models
    models <- list(dich = dich,
                   terc = terc,
                   continuous = continuous
    )
    
    
    # saveRDS( file = paste0(resultdir, '/models.rds'), models  )
    
    #check model assumptions...
    modeldiagnosticsfile <- paste0(resultdir, '/modeldiagnostics_stratum', stratum, '.pdf')
    pdf(modeldiagnosticsfile, height = 10, width = 10)
    for(model in models){
      
      print( ggcoxzph( cox.zph(model) ) )
      
    }
    dev.off()
    
    #save plots
    plots <- list(gdich, gterc)
    
    #saveRDS( file = paste0(resultdir, '/plots.rds'), plots   )
    
    ggsave( paste0(resultdir, '/gdich_stratum', stratum, '.jpg'), gdich , height = 7, width = 7)
    ggsave( paste0(resultdir, '/gterc_stratun', stratum, '.jpg'), gterc , height = 7, width = 7)
    
    #add all to lists
    results <- list(plots = plots,
                    models = models)
    
    stratumlist[[stratum]] <- results
    
    
  }
  
  
  #add stratum results to overal res
  outputlist[[scorename]] <- stratumlist
  
  
  
}




#summary

summaryfile <- paste0(suboutdir, '/summary-plots_models.pdf')
pdf(file = summaryfile, height = 9,width = 9)

for(score in names(outputlist) ){
  
  stratumlist <- outputlist[[score]]
  
  #make separate PDF page for each stratum
  for(stratum in names(stratumlist) ){
    results <- stratumlist[[stratum]]
    
    plots <- results$plots
    models <- results$models
    
    #get the coefficients, hazard ratios, confidence intervals, and P values
    summarytable <- bind_rows( lapply(models, function(x){as.data.frame( summary(x)$coefficients ) }) )
    summarytable <- bind_cols(summarytable,
                              bind_rows( lapply(models, function(x){as.data.frame( summary(x)$conf.int ) }) )[,3:4] )
    
    
    #rename the continious...
    rownames <- rownames(summarytable)
    
    #rownames[grepl('^continuous', rownames)] <- c('Univariable\nContinuous', 'Multivariable\nContinuous')
    
    rownames(summarytable) <- rownames
    
    #put conf int next to coef
    summarytable <- summarytable[,c(1,6,7,2,3,4,5)]
    
    
    #add significance marks
    summarytable$significance <- ''
    summarytable[summarytable$`Pr(>|z|)` < 0.1, "significance"] <- '.'
    summarytable[summarytable$`Pr(>|z|)` < 0.05, "significance"] <- '*'
    summarytable[summarytable$`Pr(>|z|)` < 0.01, "significance"] <- '**'
    summarytable[summarytable$`Pr(>|z|)` < 0.001, "significance"] <- '***'
    
    #round table
    summarytable[,-ncol(summarytable)] <- round( summarytable[,-ncol(summarytable)],2 )
    
    # plot it all together
    updatedplot <- (plots[[1]] + plots[[2]]) / ( gridExtra::tableGrob(summarytable)) + patchwork::plot_annotation(score) + 
      plot_layout(heights = c(1,2))
    
    
    print(updatedplot)
    
    
  }
  
  
}
dev.off()










beepr::beep()




