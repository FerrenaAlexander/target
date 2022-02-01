library(tidyverse)
library(preprocessCore)
library(survival)
library(survminer)
library(DESeq2)

set.seed(2020)

setwd("/Users/ferrenaa/Documents/deyou/target/")




r2 <- read.csv("~/Dropbox/SKP2_Data/GSE33382_R2_Phenotype.csv")


r2$Overall_Time_years <- r2$Overall_Time / 12

hist(r2$Overall_Time_years)





r2$status_code <- ifelse(r2$Deceased == 'TRUE', yes=2, no = 1)

fit <- survfit(Surv(Overall_Time_years, status_code) ~ 1, data = r2)


ggsurvplot(fit,
           pval = TRUE, conf.int = F,
           linetype = "strata", 
           #surv.median.line = "hv", 
           #risk.table = T, cumevents = T, cumcensor = T, tables.height = 0.1,
           ggtheme = theme_light())+ xlab('Time (Years)')
