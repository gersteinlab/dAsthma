### GSEA over hidden units

library(gplots)
library(RColorBrewer)
library(distances)
library(preprocessCore)
library(tsne)
library(randomForest)
library(caret)
library(pROC)
library(glmnet)
library(e1071)
library(fgsea)

load('data.RData')

### Load gene sets for GSEA
gene.set.kegg <- gmtPathways("data/gene_set/c2.cp.kegg.v6.2.symbols.gmt")
gene.set.reactome <- gmtPathways("data/gene_set/c2.cp.reactome.v6.2.symbols.gmt")
gene.set.bp <- gmtPathways("data/gene_set/c5.bp.v6.2.symbols.gmt")

###### Gene selection and GSEA
### GSEA for Hsig
# Consider the model weight over genes as the "rank"
plot.gsea <- function(gsea, pathways, ranks, path=NULL, n=10){
  if(length(ranks)==0){return(0)}
  if(!is.null(path)){pdf(path, width = 12, height = 10)}
  topPathwaysUp = gsea[ES > 0][head(order(pval), n=n), pathway]
  topPathwaysDown = gsea[ES < 0][head(order(pval), n=n), pathway]
  topPathways = c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(pathways[topPathways], ranks, gsea, 
                gseaParam = 0.5)
}

for(i in hsig){
  v.id = paste("V",i,sep="")
  wv = wt[i,]
  names(wv) = probes.gene.names
  
  pdf(paste("Figures/",v.id,".pdf", sep=""), width=6, height=4)
  plotEnrichment(gene.set.kegg[["KEGG_ASTHMA"]], wv) + labs(title=v.id)
  dev.off()
}

### Top weighted genes for each Hsig
# Remove ribosomal genes
probes.ribosomal = read.table("data/probes_ribosomal.txt")$V1
idx.non.ribosomal = (probe.id %in% probes.ribosomal == FALSE)

wt.non.ribosomal = wt[,idx.non.ribosomal]

top.gene = list()
top.probe = list()
top.gene.GO = list()

bot.gene = list()
bot.probe = list()
bot.gene.GO = list()

gsea.list = list()

N = length(probe.id)
for(i in hsig){
  v.id = paste("V",i,sep="")
  wv = wt.non.ribosomal[i,]
  names(wv) = probes.gene.names[idx.non.ribosomal]
  #wv = wt[i,var.probe.idx]
  meanv = mean(wv)
  sdv = sd(wv)
  
  #top.id = var.probe.idx[which((wv-meanv)/sdv>=2)]
  #top.id = which((wv-meanv)/sdv>=3)
  top.id = order(wv, decreasing = T)[1:200]
  top.probe[[v.id]] = probe.id[idx.non.ribosomal][top.id]
  top.gene[[v.id]] = probes.gene.names[idx.non.ribosomal][top.id]
  
  #bot.id = var.probe.idx[which(-(wv-meanv)/sdv>=2)]
  #bot.id = which(-(wv-meanv)/sdv>=3)
  bot.id = order(wv)[1:200]
  bot.probe[[v.id]] = probe.id[idx.non.ribosomal][bot.id]
  bot.gene[[v.id]] = probes.gene.names[idx.non.ribosomal][bot.id]
  
  gsea.list[[v.id]] = fgsea(gene.set.kegg, wv, nperm=1000)
  plot.gsea(gsea.list[[v.id]], gene.set.kegg, wv,paste("genes_hidden_units/GSEA/", v.id, ".pdf", sep=""))
  dev.off()
  top.gene[[v.id]] = gene.id[top.id]
  write(gene.id[top.id], file=paste("genes_hidden_units/H",i,"_top.txt",sep=""))
  write(gene.id[bot.id], file=paste("genes_hidden_units/H",i,"_bot.txt",sep=""))
}
