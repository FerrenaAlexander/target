library(tidyverse)




### figs for D Geller grant 2022.05.11 #####
# survival analysis of macs in target --> done check other follders
# bulk RNAseq: TME / macs big partof TME --> will show in target, check in kuijjer
# TKO vs DKO big change in macs infiltratio --> shown in paper figures
# scrnaseq, macs... need to think about this...
# - for this, repeat mac specturm analysis....




# target show macs and TME

# links to a live github: 2022.01.05
mcp <- readRDS("data/mcpresults_2022.rds")

#exclude fibs, super high, maybe confounded by sarcoma...
mcp <- mcp[rownames(mcp)!='Fibroblasts',]

#melt it for plotting
x <- reshape2::melt(mcp)

#plot
mcp_plot <- ggplot(x, aes(Var1, value))+geom_boxplot()+ 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab('Cell type') + ylab('Microevironment Cell Popullation Counter score')




# using cibersort...
cbx <- read.csv('data/cibersortresults/CIBERSORTx_Job2_Results.csv')

#keep same
sampname <- str_split_fixed(cbx$Mixture, '-', Inf)[,1:3]
sampname <- apply(sampname, MARGIN = 1, FUN = function(x){paste(x, collapse ='-')})
cbx$sample <- sampname
cbx <- cbx[cbx$sample %in% colnames(mcp),]
cbx <- cbx[match(colnames(mcp), cbx$sample),]

#cor them...
mcpt <- as.data.frame(t(mcp))

cor.test(cbx$Monocyte.Macrophage, mcpt$`Monocytic lineage`)

#plot them...
cbx <- cbx[,!(colnames(cbx) %in% c('Mixture', 'Tumor', "MSC.Fibroblast" ,  "Plasma", "P.value", "Correlation", "RMSE", "Absolute.score..sig.score.", 'sample' ))]

y <- reshape2::melt(cbx)


cbx_plot <- ggplot(y, aes(variable, value ))+geom_boxplot()+ 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab('Cell type') + ylab('CIBERSORTx score, trained with DKO scRNA-seq')

mcp_plot
cbx_plot




comb <- data.frame(MCP_monoctye = mcpt$`Monocytic lineage`,
                   cbx_monocyte = cbx$Monocyte.Macrophage,
                   cbx_osteocllast = cbx$Osteoclast
                   
                   )


dir.create('geller_grant')

pdf('geller_grant/monocytes_percent_TME_target.pdf', height = 5, width = 5)

mcp_plot
cbx_plot

dev.off()



