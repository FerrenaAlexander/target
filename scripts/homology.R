library(tidyverse)
library(preprocessCore)
library(survival)
library(survminer)
library(DESeq2)
library(patchwork)

library(biomaRt) #mouse-human homologs

library(MCPcounter) #for calculating microenv scores in human

set.seed(2020)



### final followup from tko vs dko bulk rnaseq paper

#Multivar survival model 
# with both monocyte MCP score AND 
# TKO overexp gene list


#Subtract myeloid genes from TKO: 
# remove DE genes present in msigdb myeloid pathways, 
# check survival








#### prep run #####


#set up biomart function for mouse to human
human = NULL; mouse = NULL

while(is(human) != 'Mart'){ human = try( useMart("ensembl", dataset = "hsapiens_gene_ensembl",  host = "http://feb2021.archive.ensembl.org/")) }

while(is(mouse) != 'Mart'){ mouse = try( useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "http://feb2021.archive.ensembl.org/")) }




### try to get all human gene mappings...

gencode <-read.csv('~/Dropbox//data/bangdata/2021october-TKOvsDKO/data/gencode.vM23.ids_names_types.csv')

gencode$id_with_version <- gencode$gene_id
split <- str_split_fixed(string = gencode$gene_id, '\\.', 2)
gencode$gene_id <- split[,1] ; rm(split)


genes <- gencode$gene_id
genesV2 = getLDS(attributes = c('ensembl_gene_id',"mgi_symbol"), 
                 filters = "ensembl_gene_id", 
                 values = genes , 
                 mart = mouse, 
                 attributesL = c('ensembl_gene_id',"hgnc_symbol"), 
                 martL = human, 
                 uniqueRows=T)



#try to eliminate all dups....
nodups <- genesV2


first <- names(table(nodups[,1])[table(nodups[,1])==1])
second <- names(table(nodups[,2])[table(nodups[,2])==1])
third <- names(table(nodups[,3])[table(nodups[,3])==1])
fourth <- names(table(nodups[,4])[table(nodups[,4])==1])

nodups <- nodups[nodups[,1] %in% first,]
nodups <- nodups[nodups[,2] %in% second,]
nodups <- nodups[nodups[,3] %in% third,]
nodups <- nodups[nodups[,4] %in% fourth,]


write.csv('data/biomart_nodups_may05-2022_feb2021archive.csv', x = nodups, 
          row.names = F, quote = F)
