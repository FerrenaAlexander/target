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





#age
clin$Age_At_Diagnosis_Years <- clin$`Age at Diagnosis in Days` / 365

ggplot(clin, aes(Age_At_Diagnosis_Years))+
  geom_histogram(col='black', fill='steelblue')+
  theme_classic()



#met
clin$Primary_Metastasis <- "Non-metastatic"
clin[!grepl(pattern = 'Non', x = clin$`Disease at diagnosis`),"Primary_Metastasis"] <- 'Metastatic'


ggplot(clin, aes(fill = Primary_Metastasis, x = Primary_Metastasis))+
  geom_bar()+
  theme_linedraw()+
  scale_fill_brewer(palette = 'Set1')



#sex
ggplot(clin, aes(fill = Gender, x = Gender))+
  geom_bar()+
  theme_linedraw()
