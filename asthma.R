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


dat.raw <- read.table("data/sputum_combat_RIN_adjusted_all.txt",sep="\t")
dat <- normalize.quantiles(as.matrix(dat.raw))
colnames(dat) <- colnames(dat.raw)

frq <- read.table("data/freq.sample.clinic.txt")
non.ctrl=which(frq$V9!="CONTROL")
sev <- factor(frq$V9[non.ctrl], labels = c("MILD", "MODERATE", "SEVERE"))

hid <- as.matrix(read.table("data/New.sput.cycle100.h50.lr0.1.err0.01.hidden.txt",sep="\t"))

wt.raw <- as.matrix(read.table("data/New.sput.cycle100.h50.lr0.1.err0.01.weight.txt",sep="\t"))
wt.raw <- wt.raw[,-1]
probe.id <- rownames(read.table("data/sputum_combat_RIN_adjusted_all.txt",sep="\t"))
colnames(wt.raw) <- probe.id

rownames(dat) <- probe.id

wt <- wt.raw
#wt <- t(normalize.quantiles(t(wt.raw)))
#rownames(wt) <- rownames(wt.raw)
#colnames(wt) <- colnames(wt.raw)

tea <- read.table("data/TEA_kegg.dist_kmeans_3.txt")
tea <- as.factor(paste("TEA",tea$V2, sep=""))

hid0 = hid[non.ctrl,-1]
colnames(hid0) <- paste("H", seq(1,50), sep="")

#colnames(hid0) = paste("V",seq(1,50),sep="")
hvar = which(apply(hid0, 2, var)>0.001) # id of hidden units with variance

tmp.pal=brewer.pal(5,"Set2");
names(tmp.pal)=c("TEA1","TEA2","TEA3","ASTHMA");

tmp2.pal=tmp.pal
names(tmp2.pal)<-c("CONTROL","MILD","MODERATE","SEVERE");

c1 <- tmp.pal[as.character(tea)]
c2 <- tmp2.pal[as.character(sev)];

### Heatmap
pdf("hvar.pdf")
par(xpd=TRUE)
heatmap.2(hid0[,hvar], trace="none", cexCol=1, col=colorRampPalette(c("white","red"))(32), scale="none",dendro="row",srtCol=270,offsetCol=0,offsetRow=0,adjCol=c(0.5,1),cexRow=0.5,margin=c(3,3),RowSideColors=c1,density.info="none",labRow=NA,key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Hidden variables")
legend(x=0.9,y=1.1,c("TEA1","TEA2","TEA3"),title="Side Color",fill=tmp.pal[c("TEA1","TEA2","TEA3")],cex=0.6)
dev.off()

### Affinity matrix
dm <- as.matrix(distances(hid0[,hvar]))
am <- exp(-dm^2/2)
pdf("Figures/hvar_affinity.pdf")
par(xpd=TRUE)
heatmap.2(am, trace="none", cexCol=0.5, col=colorRampPalette(c("white","red"))(32), scale="none",dendro="row",srtCol=270,offsetCol=0,offsetRow=0,cexRow=0.5,margin=c(3,3),RowSideColors=c1,ColSideColors=c1,density.info="none",key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Affinity")
legend(x=0.9,y=1.1,c("TEA1","TEA2","TEA3"),title="Side Color",fill=tmp.pal[c("TEA1","TEA2","TEA3")],cex=0.6)
dev.off()

### Correlation with cluster
n = length(hid0[1,])

pdf("Figures/spcorr_all.pdf", width=10, height=6)
par(mfrow=c(4,6))
par(mar=c(3,1,3,1))

spcor = rep(0, n)
    
for(i in hvar){
  boxplot(hid0[,i]~tea, main=paste("V",i,sep=""), cex=0.5, ylim=c(-0.2,1),col=tmp.pal[1:3])
  spcor[i] = cor(hid0[,i],as.numeric(tea), method="spearman")
  text(x=2, y=-0.11, paste("Spearman correlation = ",round(spcor[i],3), sep=""),cex=0.8)
  
}

dev.off()

hsig = which(abs(spcor)>0.65) # id of hidden units with significant correlation with tea

pdf("Figures/spcorr.pdf", width=10, height=6)
par(mfrow=c(2,3))
par(mar=c(3,1,3,1))

for(i in hsig){
  if(i==36){plot.new()}
  boxplot(hid0[,i]~tea, main=paste("H",i,sep=""), cex=0.5, ylim=c(-0.2,1),col=tmp.pal[1:3])
  text(x=2, y=-0.11, paste("Spearman correlation = ",round(spcor[i],3), sep=""),cex=1)
}
dev.off()

### Heatmap
pdf("Figures/hsig.pdf")
par(xpd=TRUE)
heatmap.2(hid0[,hsig], trace="none", cexCol=1, col=colorRampPalette(c("white","red"))(32), scale="none",dendro="row",srtCol=270,offsetCol=0,offsetRow=0,adjCol=c(0.5,1),cexRow=0.5,margin=c(3,3),RowSideColors=c1,density.info="none",labRow=NA,key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Hidden variables")
legend(x=0.9,y=1.1,c("TEA1","TEA2","TEA3"),title="Side Color",fill=tmp.pal[c("TEA1","TEA2","TEA3")],cex=0.6)
dev.off()

### Affinity matrix
pdf("Figures/hsig_affinity.pdf")
dm <- as.matrix(distances(hid0[,hsig]))
am <- exp(-dm^2/2)
par(xpd=TRUE)
heatmap.2(am, trace="none", cexCol=0.5, col=colorRampPalette(c("white","red"))(32), scale="none",dendro="row",srtCol=270,offsetCol=0,offsetRow=0,cexRow=0.5,margin=c(3,3),RowSideColors=c1,ColSideColors=c1,density.info="none",key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Affinity")
legend(x=0.9,y=1.1,c("TEA1","TEA2","TEA3"),title="Side Color",fill=tmp.pal[c("TEA1","TEA2","TEA3")],cex=0.6)
dev.off()

### Weight distribution
pdf("Figures/weights.pdf", width=10, height=6)
par(mfrow=c(2,3))
par(mar=c(3,1,3,1))

for(i in hsig){
  if(i==36){plot.new()}
  wv = wt[i,]
  meanv = mean(wv)
  sdv = sd(wv)
  plot(density(wv), main=paste("H",i,sep=""), cex=0.5)
  lines(rep(meanv-sdv,2), c(0,max(density(wv)$y)), lty=3, col="blue")
  lines(rep(meanv+sdv,2), c(0,max(density(wv)$y)), lty=3, col="blue")
}
dev.off()

h.col = c("red","red","deepskyblue","deepskyblue","deepskyblue")
h.lty = c(2,1,4,2,1)

pdf("Figures/weights_combined.pdf", width=8, height=6)
plot(1, xlim=c(-0.13,0.5), ylim=c(0,40),xlab='',ylab='density')
for(i in seq(1,length(hsig))){
  wv = wt[hsig[i],]
  lines(density(wv), col = h.col[i], lty=h.lty[i])
}
legend("topright",paste("H",hsig,sep=""),col=h.col, lty=h.lty)
dev.off()

### Load annotation
myinf2 = "data/GPL6244.annot"

## bread in the gene name for each probeset
conIn = file(myinf2, "r")
data = readLines(conIn, -1)
close(conIn)
mysta = grep("platform_table_begin", data)
myend = grep("platform_table_end", data)
data = read.table(myinf2, sep="\t", header=T, row.names=1, skip=mysta, nrows=myend-mysta-2, quote="", comment.char ="")
tmp = row.names(data)

#GPL = data[, "Gene.symbol"]
GPL = data[, "GenBank.Accession"]
GS = data[, "Gene.symbol"]
GT = data[,"Gene.title"]
names(GPL) = tmp
names(GS) = tmp
names(GT) = tmp

gene.id = as.character(GPL[probe.id])



################## Gene selection and GSEA
sd.dat = apply(dat, 1, sd) # variance of expression
sd.wt = apply(wt, 2, sd) # variance of expression
par(mfrow=c(2,3))

for(i in hsig){
  if(i==36){plot.new()}
  wv = wt[i,]
  meanv = mean(wv)
  sdv = sd(wv)
  plot(wv~apply(dat.raw, 1, mean), main=paste("H",i,sep=""), cex=0.5, xlab='mean expression',ylab='weight')
  #lines(rep(meanv-sdv,2), c(0,max(density(wv)$y)), lty=3, col="blue")
  #lines(rep(meanv+sdv,2), c(0,max(density(wv)$y)), lty=3, col="blue")
}

probes.gene.names = unlist(lapply(as.character(GS[probe.id]), function(x){strsplit(x, "/")[[1]][1]}))
names(probes.gene.names) = probe.id
gene.set.kegg <- gmtPathways("data/gene_set/c2.cp.kegg.v6.2.symbols.gmt")
gene.set.reactome <- gmtPathways("data/gene_set/c2.cp.reactome.v6.2.symbols.gmt")
gene.set.bp <- gmtPathways("data/gene_set/c5.bp.v6.2.symbols.gmt")

### Top scored genes
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


plot.gsea <- function(gsea, pathways, ranks, path=NULL, n=10){
  if(length(ranks)==0){return(0)}
  if(!is.null(path)){pdf(path, width = 12, height = 10)}
  topPathwaysUp <- gsea[ES > 0][head(order(pval), n=n), pathway]
  topPathwaysDown <- gsea[ES < 0][head(order(pval), n=n), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(pathways[topPathways], ranks, gsea, 
                gseaParam = 0.5)
}

### GSEA
for(i in hsig){
  v.id = paste("V",i,sep="")
  wv = wt[i,]
  names(wv) = probes.gene.names
  
  pdf(paste("Figures/",v.id,".pdf", sep=""), width=6, height=4)
  plotEnrichment(gene.set.kegg[["KEGG_ASTHMA"]], wv) + labs(title=v.id)
  dev.off()
}

### Top genes
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

top.probe.list = unique(unlist(top.probe))
bot.probe.list = unique(unlist(bot.probe))
probe.list = unique(c(top.probe.list,bot.probe.list))
#probe.list = top.probe.list

#x = t(dat[rel.genes,non.ctrl])
x = t(dat[top.probe.list,non.ctrl])

### Heatmap of selected gene weight
# set labels
neg.top.list = unique(c(top.probe$V26,top.probe$V27))
pos.top.list = unique(c(top.probe$V36,top.probe$V38, top.probe$V45))

shared.list = intersect(neg.top.list,pos.top.list)
neg.unique.list = setdiff(neg.top.list,intersect(neg.top.list,pos.top.list))
pos.unique.list = setdiff(pos.top.list,intersect(neg.top.list,pos.top.list))

write(as.character(GPL[pos.unique.list]), file=paste("genes_hidden_units/top.pos_unique.txt",sep=""))
write(as.character(GPL[neg.unique.list]), file=paste("genes_hidden_units/top.neg_unique.txt",sep=""))


probe.set.labels = rep("", length(top.probe.list))
names(probe.set.labels) = c(shared.list,neg.unique.list,pos.unique.list)

probe.set.labels[shared.list] = "shared"
probe.set.labels[neg.unique.list] = "neg-unique"
probe.set.labels[pos.unique.list] = "pos-unique"

set.pal=c("green","blue","orange")
names(set.pal)=c("shared","neg-unique","pos-unique");
set.c = set.pal[probe.set.labels]

sev.pal=brewer.pal(5,"Set2");
names(sev.pal)<-c("CONTROL","MILD","MODERATE","SEVERE");

sev.c <- sev.pal[as.character(sev_sub)];

par(xpd=TRUE)
heatmap.2(t(wt[hsig,c(shared.list,neg.unique.list,pos.unique.list)]), trace="none", cexCol=1, 
          col=colorRampPalette(c("white","white","red","red"))(256), 
          RowSideColors=set.c, 
          labCol = hsig, Rowv="none",Colv="none",dendrogram="none")
legend(x=0.8,y=1.1,c("Shared","Negative set-unique","Positive set-unique"),title="Side Color",fill=set.pal[c("shared","neg-unique","pos-unique")],cex=0.6)

par(xpd=TRUE)
heatmap.2(as.matrix(dat[c(shared.list,neg.unique.list,pos.unique.list),sev_sub_id]), trace="none", cexCol=1, 
          col=colorRampPalette(c("blue","white","red"))(256), 
          RowSideColors=set.c,  
          Rowv="none",dendrogram="none", scale="column")
legend(x=0.8,y=1.1,c("Shared","Negative set-unique","Positive set-unique"),title="Side Color",fill=set.pal[c("shared","neg-unique","pos-unique")],cex=0.6)


### Gene expression analysis
tea1_deg_probe <- as.character(read.table("data/tea_degs/tea1_degs_id.txt")$V1)
tea3_deg_probe <- as.character(read.table("data/tea_degs/tea3_degs_id.txt")$V1)
deg = unique(c(tea1_deg_probe, tea3_deg_probe))

### Intersection between gene lists
intersect_mat = matrix(0L, nrow = length(hsig), ncol = length(hsig))

for(i in seq(1,length(hsig))){
  for(j in seq(1,length(hsig))){
    intersect_mat[i,j] = length(intersect(top.probe[[i]],top.probe[[j]]))
  }
}

write(as.character(GPL[top.probe.list]), file=paste("genes_hidden_units/top.txt",sep=""))
write(as.character(GS[top.probe.list]), file=paste("genes_hidden_units/top_symbol.txt",sep=""))
write(as.character(GPL[bot.probe.list]), file=paste("genes_hidden_units/bot.txt",sep=""))


### Heatmap of selected gene expression
par(xpd=TRUE)
heatmap.2(as.matrix(dat[probe.list,non.ctrl]), trace="none", col=colorRampPalette(c("green","black","red"))(256), scale="row",offsetCol=0,offsetRow=0,adjCol=c(0.5,1),margin=c(3,3),ColSideColors=c1,density.info="none",labRow=NA, labCol=NA,key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Expression")
legend(x=0.9,y=1.1,c("TEA1","TEA2","TEA3"),title="Side Color",fill=tmp.pal[c("TEA1","TEA2","TEA3")],cex=0.6)

### Correlation of weights of top-weighted genes
wt.cor = matrix(0L, nrow = length(hsig), ncol = length(hsig))

for(i in seq(1,length(hsig))){
  for(j in seq(i,length(hsig))){
    wt.cor[i,j] = cor(wt[i,top.probe.list], wt[j,top.probe.list],method = "pearson")
    wt.cor[j,i] = wt.cor[i,j]
  }
}

### TSNE
genes with relevant functional terms
rel.genes <- as.character(read.table("genes_hidden_units/top_genes_merged_probe.txt")$V1)
rel.genes <- rel.genes[rel.genes %in% rownames(dat)]

for(p in c(5,15,25,35,45)){
  x_tsne = tsne(x[sev_sub_id,imp.probes],perplexity=p)
  plot(x_tsne[,1], x_tsne[,2], col=sev_sub, xlab="TSNE1",ylab="TSNE2")
}

### PCA
x_pca = prcomp(t(x[,imp.probes]))$rotation
plot(x_pca[sev_sub_id,1], x_pca[sev_sub_id,2], col=sev_sub, xlab="PC1",ylab="PC2")

par(xpd=TRUE)
heatmap.2(t(x), trace="none", col=colorRampPalette(c("green","black","red"))(256), scale="row",offsetCol=0,offsetRow=0,adjCol=c(0.5,1),margin=c(3,3),ColSideColors=c1,density.info="none",labRow=NA, labCol=NA,key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Expression")
legend(x=0.9,y=1.1,c("TEA1","TEA2","TEA3"),title="Side Color",fill=tmp.pal[c("TEA1","TEA2","TEA3")],cex=0.6)

################## Prediction of clusters
### Random forest
rf_performance <- function(x, y){
  id.neg = which(as.numeric(y)==1)
  id.pos = which(as.numeric(y)==2)
  N.neg = sum(as.numeric(y)==1)
  N.pos = sum(as.numeric(y)==2)
  acc_all = rep(0,10)
  dat.all = data.frame(label=y,x)
  
  for(n in seq(1,10)){
    predicted.probs = NULL
    true.lables = NULL
    id.parts.neg = sample(1:4,size=N.neg,replace=TRUE)
    id.parts.pos = sample(1:4,size=N.pos,replace=TRUE)
    acc = rep(0,4)
    for(i in seq(1,4)){
      id = c(id.neg[which(id.parts.neg==i)],id.pos[which(id.parts.pos==i)])
      dat.train = dat.all[-id,]
      dat.test = dat.all[id,]
      rf = train(label~., dat.train)
      predicted.prob = predict(rf, dat.test, type="prob")[,"TRUE"]
      
      predicted.probs = c(predicted.probs, predicted.prob)
      true.lables = c(true.lables, dat.test$label)
      
      #result.roc = roc(dat.test$label, predicted.prob[,"TRUE"]) # Draw ROC curve.
      #cm = confusionMatrix(y[id], predict(rf, x[id,]))
      #acc[i] = cm$overall[1]
      #acc[i] = pROC::auc(result.roc)
    }
    result.roc = roc(true.lables, predicted.probs) # Draw ROC curve.
    acc_all[n] = pROC::auc(result.roc)
    #acc_all[n] = mean(acc)
  }
  return(mean(acc_all))
#    return(mean(acc))
}

plot_rf_performance <- function(x, y, path=NULL, main=NULL, add=F, col="red"){
  if(!is.null(path)){pdf(path)}
  id.neg = which(as.numeric(y)==1)
  id.pos = which(as.numeric(y)==2)
  N.neg = sum(as.numeric(y)==1)
  N.pos = sum(as.numeric(y)==2)
  dat.all = data.frame(label=y,x)
  
  predicted.probs = NULL
  true.lables = NULL
  id.parts.neg = sample(1:4,size=N.neg,replace=TRUE)
  id.parts.pos = sample(1:4,size=N.pos,replace=TRUE)

  for(i in seq(1,4)){
    id = c(id.neg[which(id.parts.neg==i)],id.pos[which(id.parts.pos==i)])
    dat.train = dat.all[-id,]
    dat.test = dat.all[id,]
    rf = train(label~., dat.train)
    predicted.prob = predict(rf, dat.test, type="prob")[,"TRUE"]
      
    predicted.probs = c(predicted.probs, predicted.prob)
    true.lables = c(true.lables, dat.test$label)
  }
  result.roc = roc(true.lables, predicted.probs) # Draw ROC curve.
  plot(result.roc, col=col,main=main, add=add, print.auc=T)
}


tea_1and3_id = c(which(tea=="TEA1"), which(tea=="TEA3"))
sev_sub_id = c(which(sev=="MILD"), which(sev=="SEVERE"))
sev_sub = as.factor(as.character(sev[sev_sub_id]))

t.tea = as.factor(as.numeric(tea[tea_1and3_id])==3)
t.sev =  as.factor(as.numeric(sev[sev_sub_id])==3)

h.sig = hid0[, hsig]
h.var = hid0[, hvar]

rf_performance(h.sig[tea_1and3_id, ], t.tea)
rf_performance(h.var[tea_1and3_id, ], t.tea)
rf_performance(x[tea_1and3_id, ], t.tea)

rf_performance(h.sig[sev_sub_id, ], t.sev)
rf_performance(h.var[sev_sub_id, ], t.sev)
rf_performance(x[sev_sub_id, ], t.sev)

################## Importance of genes
gene.imp = data.frame(n=rep(0,length(x[1,])))
rownames(gene.imp) = colnames(x)

for(i in seq(1,50)){
  dat.train = data.frame(label=t.sev,x=x[sev_sub_id, ])
  colnames(dat.train)[-1] = colnames(x)
  rf = train(label~., dat.train)
  var.imp = varImp(rf)$importance
  gene.imp = cbind(gene.imp, as.numeric(var.imp$Overall))
}

gene.imp = gene.imp[,-1]
colnames(gene.imp) = paste("rf",seq(1,50), sep="")
gene.imp.mean = rowMeans(gene.imp)
names(gene.imp.mean) = probes.gene.names[names(gene.imp.mean)]

gsea = fgsea(gene.set.kegg, gene.imp.mean, nperm=1000)
plot.gsea(gsea, gene.set.kegg, gene.imp.mean)

N.feat = length(x[1,])
mean.imp = mean(gene.imp.mean)
std.imp = sd(gene.imp.mean)
#imp.probes = names(gene.imp.mean)[(gene.imp.mean-mean.imp)/std.imp > 1]
imp.probes = rownames(gene.imp)[order(gene.imp.mean)[(N.feat-49):N.feat]]
top.imp.probes = rownames(gene.imp)[order(gene.imp.mean)[(N.feat-9):N.feat]]

heatmap.2(as.matrix(dat[top.imp.probes,sev_sub_id]), trace="none", col=colorRampPalette(c("green","black","red"))(256), scale="row",offsetCol=0,offsetRow=0,adjCol=c(0.5,1),margin=c(3,3),ColSideColors=c3,density.info="none",labRow=NA, labCol=NA,key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Expression")

write(as.character(GPL[imp.probes]), file=paste("clinical/sev_important_genes.txt",sep=""))
write(as.character(GS[imp.probes]), file=paste("clinical/sev_important_genes_symbol.txt",sep=""))

rf_performance(x[sev_sub_id, imp.probes], t.sev)

# Compare with random
rf_rand = rep(0,10)
for(i in seq(1,10)){
  rand.probes = sample(colnames(x),50)
  rf_rand[i] = rf_performance(x[sev_sub_id, rand.probes], t.sev)
}
mean(rf_rand)


### Pval of expression
probes.pval = rep(0, length(top.probe.list))

for(i in seq(1,length(top.probe.list))){
  v = dat[top.probe.list[i],sev_sub_id]
  
  probes.pval[i] = t.test(v[sev_sub=="MILD"], v[sev_sub=="SEVERE"], alternative = "two.sided")$p.value
}

### Plot of importance
plot(density(gene.imp.mean))

pdf("Figures/top.probes.importance.pdf", width=10, height=4)
imp.probes.id  = order(gene.imp.mean, decreasing = T)[1:50]
gene.imp.std = apply(gene.imp, 1, FUN=sd)

bp = barplot(gene.imp.mean[imp.probes.id], names.arg=NA,las=2, ylim=c(0,104), col="blue")
arrows(bp, gene.imp.mean[imp.probes.id]-gene.imp.std[imp.probes.id], 
       bp, gene.imp.mean[imp.probes.id]+gene.imp.std[imp.probes.id], 
       length=0.02, angle=90, code=3)
text(bp, -0.5, srt = 60, adj= 1, xpd = TRUE, labels = GS[imp.probes], cex=0.8)
dev.off()

### Expression and relationship with severity
tmp3.pal=tmp.pal
names(tmp3.pal)<-c("CONTROL","MILD","MODERATE","SEVERE");
c3 <- tmp3.pal[sev_sub]

pdf("Figures/top.probes.expression.pdf", width=10, height=6)
par(mfrow=c(3,4))
par(mar=c(3,1,3,1))

for(i in top.imp.probes){
  v = dat[i,sev_sub_id]
  boxplot(v~sev_sub, main=GS[i], cex=0.5,col=tmp.pal[1:3])
  t.test.pval = t.test(v[sev_sub=="MILD"], v[sev_sub=="SEVERE"], alternative = "two.sided")$p.value
  text(x=1.5, y=(max(v) - min(v))*0.1+min(v), paste("T-test pval =\n",round(t.test.pval,3), sep=""),cex=1)
}

dev.off()

### Genes related to autoimmune response, etc.
pdf("Figures/probes.immune.expression.pdf", width=10, height=6)
par(mfrow=c(2,3))
par(mar=c(3,1,3,1))

for(i in c("7905571","8117800","8156228","8178891","8177732")){
  v = dat[i,sev_sub_id]
  boxplot(v~sev_sub, main=GS[i], cex=0.5,col=tmp.pal[1:3])
  t.test.pval = t.test(v[sev_sub=="MILD"], v[sev_sub=="SEVERE"], alternative = "two.sided")$p.value
  text(x=1.5, y=(max(v) - min(v))*0.1+min(v), paste("T-test pval =\n",round(t.test.pval,3), sep=""),cex=1)
}

dev.off()


#dat.train = data.frame(label=t.sev,x=x[sev_sub_id,imp.probes])
#colnames(dat.train)[-1] = imp.probes
#rf = train(label~., dat.train)

### Correlation with severity
n = length(imp.probes)

par(mfrow=c(4,5))
par(mar=c(2,2,2,2))

wilcox.pval = rep(1, n)

for(i in imp.probes){
  expi = x[sev_sub_id, i]
  boxplot(expi~sev_sub, main=GS[i], cex=0.5, ylim=c(min(expi), max(expi)))
  wilcox.pval[i] = wilcox.test(expi[sev_sub=="MILD"], expi[sev_sub=="SEVERE"], exact = F, alternative = "greater")$p.value
  text(x=1.5, y=min(expi)+0.1*(max(expi)-min(expi)), paste("Wilcox pval =\n",wilcox.pval[i]),cex=0.8)
}


### Correlation with severity

pdf("Figures/sev_spcorr.pdf", width=10, height=6)
par(mfrow=c(2,3))
par(mar=c(3,1,3,1))

for(i in hsig){
  if(i==36){plot.new()}
  h = hid0[,i]
  spcor0 = cor(hid0[,i],as.numeric(sev), method="spearman")
  boxplot(hid0[,i]~sev, main=paste("H",i,sep=""), cex=0.5, ylim=c(-0.2,1),col=tmp.pal[1:3])
  text(x=2, y=-0.11, paste("Spearman correlation = ",round(spcor0,3), sep=""),cex=1)
}
dev.off()

n = length(hid0[1,])

pdf("Figures/sev_spcorr_all.pdf", width=10, height=6)
par(mfrow=c(4,6))
par(mar=c(3,1,3,1))

for(i in hvar){
  h = hid0[,i]
  spcor0 = cor(hid0[,i],as.numeric(sev), method="spearman")
  boxplot(hid0[,i]~sev, main=paste("H",i,sep=""), cex=0.5, ylim=c(-0.2,1),col=tmp.pal[1:3])
  text(x=2, y=-0.11, paste("Spearman correlation = ",round(spcor0,3), sep=""),cex=0.8)
}

dev.off()


### Plot ROC
plot_rf_performance(t(dat[deg, sev_sub_id]), t.sev, path="roc.pdf", col="blue")
plot_rf_performance(h.sig[sev_sub_id, ], t.sev, col="purple", add=T)
plot_rf_performance(h.var[sev_sub_id, ], t.sev, col="green", add=T)
plot_rf_performance(x[sev_sub_id, ], t.sev, col="cyan", add=T)
plot_rf_performance(x[sev_sub_id, imp.probes], t.sev, col="red", add=T)
dev.off()

