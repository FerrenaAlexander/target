library(tidyverse)
library(preprocessCore)
library(survival)
library(survminer)
library(DESeq2)
library(patchwork)

library(foreach)
library(doParallel)

library(biomaRt) #mouse-human homologs

library(MCPcounter) #for calculating microenv scores in human

set.seed(2020)





### this is to test the TKO output in the kuijjer dataset.






#set up biomart function for mouse to human
human = NULL; mouse = NULL

while(is(human) != 'Mart'){ human = try( useMart("ensembl", dataset = "hsapiens_gene_ensembl")) }

while(is(mouse) != 'Mart'){ mouse = try( useMart("ensembl", dataset = "mmusculus_gene_ensembl")) }





### mouse to human

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




### fix long names function
fixlongnames <- function(pnames){
  
  pnewnames <- c()
  for(p in pnames){
    
    
    if(str_length(p) > 50){
      
      #try to find and replace underscore closest to halfway...
      halfway <- round(str_length(p)/2)
      underscorepos <- str_locate_all(p,'_')[[1]][,1]
      
      distances <- abs(halfway-underscorepos)
      
      which_und_is_closest <- which.min(distances)
      
      split <- str_split_fixed(p,'_',Inf)
      
      newp <- paste0(paste(split[1,1:which_und_is_closest], collapse = '_'),
                     '        ', '\n',
                     paste(split[1,(which_und_is_closest+1):ncol(split)], collapse = '_')
      )
      
      
      
    } else{newp <- p}
    
    pnewnames <- c(pnewnames, newp)
    
  }
  
  pnewnames
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









#readin gem, make first col row names
gem <- read.csv('data/kuijjer/parsed-gem.csv')
rownames(gem) <- gem[,1] ; gem <- gem[,-1]



#read in kuijjer md
md <- read.csv('data/kuijjer/parsed-md.csv')



#read in genes
genes <- read.csv('data/kuijjer/parsed-genesymbols.csv')



#get only samples with survival
md <- md[!is.na(md$metastasis),]
gem <- gem[,match(md$sample, colnames(gem))]




#select genes

gene = 'ILMN_1665538'


genevec <- t(gem[gene,])[,1]

#clindf
clindf <- data.frame(sample = md$sample,
                     met = md$metastasis,
                     deceased = md$deceased,
                     tof = md$tof,
                     tom1 = md$tom1,
                     tod = md$tod,
                     gene = genevec
)



#dichotomize
clindf$genedich <- 'Lo'
clindf[clindf$gene >= median(clindf$gene), "genedich"] <- 'Hi'



#try to make combined met/censoring var
# for met=true, use the tom1
# for met = false, use tof

#one person not deceased a nd no met is missing tim, just mark is as followup=0...
md[is.na(md$tof),'tof'] <- 0

#set time to met = tof, this is for met=false
clindf$time_to_met <- md$tof

#for met=true, set time_to_met=tom1
clindf[clindf$met=='true',"time_to_met"] <- clindf[clindf$met=='true',"tom1"]


clindf$met_recode <- 'Censored'
clindf[clindf$met=='true', "met_recode"] <- 'Event'
clindf$met_recode <- factor(clindf$met_recode, levels = c('Event', 'Censored'))




df <- data.frame(status = clindf$met_recode,
                 surv = clindf$time_to_met,
                 dichotomized = clindf$genedich)
df$status_code <- ifelse(df$status == 'Censored', yes=0, no = 1)

fit <- survfit(Surv(surv, status_code) ~ dichotomized, data = df)

ggsurvplot(fit,
           conf.int = F,
           pval = TRUE, #pval.coord = c(0, 0.15),
           pval.method = T, #pval.method.coord = c(0, 0.25),
           linetype = "strata",  
           test.for.trend = F,
           ggtheme = theme_classic2(), 
           palette = c("#377EB8" , "#E41A1C") )+ 
  xlab('Time (Months)')




