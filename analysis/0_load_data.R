### Load raw expression data and hidden layer values
### Hierarchical clustering and correlation analysis

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

###### Load data
### Raw expression data & quantile normalize
dat.raw = read.table("data/sputum_combat_RIN_adjusted_all.txt",sep="\t")
dat = normalize.quantiles(as.matrix(dat.raw))
colnames(dat) = colnames(dat.raw)

### Clinical information
frq = read.table("data/freq.sample.clinic.txt")
non.ctrl=which(frq$V9!="CONTROL")
sev = factor(frq$V9[non.ctrl], labels = c("MILD", "MODERATE", "SEVERE"))

tea = read.table("data/TEA_kegg.dist_kmeans_3.txt")
tea = as.factor(paste("TEA",tea$V2, sep=""))

### Other clinical traits
clinical = read.table("data/clinical.csv",sep=",")
colnames(clinical) = c("Array_ID", "PRED.FEV1", "PRED.FVC", "POST.FEV1", "POST.FVC", "POST.FEV1/FVC", "PRE.FEV1", "PRE.FVC", "PRE.FEV1/FVC")
rownames(clinical) = clinical[,1]
clinical = clinical[,-1]

### Autoencoder hidden layer values and model weights
hid = as.matrix(read.table("data/New.sput.cycle100.h50.lr0.1.err0.01.hidden.txt",sep="\t"))

wt.raw = as.matrix(read.table("data/New.sput.cycle100.h50.lr0.1.err0.01.weight.txt",sep="\t"))
wt.raw = wt.raw[,-1]
probe.id = rownames(read.table("data/sputum_combat_RIN_adjusted_all.txt",sep="\t"))
colnames(wt.raw) = probe.id

rownames(dat) = probe.id

wt = wt.raw

hid0 = hid[non.ctrl,-1]
colnames(hid0) = paste("H", seq(1,50), sep="")
hvar = which(apply(hid0, 2, var)>0.001) # id of hidden units with higher variance


###### Visualization and clustering
### Heatmap and hierarchical clustering of hidden layer values (non-control samples only)
tmp.pal=brewer.pal(5,"Set2");
names(tmp.pal)=c("TEA1","TEA2","TEA3","ASTHMA");

tmp2.pal=tmp.pal
names(tmp2.pal)<-c("CONTROL","MILD","MODERATE","SEVERE");

c1 = tmp.pal[as.character(tea)]
c2 = tmp2.pal[as.character(sev)];

pdf("hvar.pdf")
par(xpd=TRUE)
heatmap.2(hid0[,hvar], trace="none", cexCol=1, col=colorRampPalette(c("white","red"))(32), scale="none",dendro="row",srtCol=270,offsetCol=0,offsetRow=0,adjCol=c(0.5,1),cexRow=0.5,margin=c(3,3),RowSideColors=c1,density.info="none",labRow=NA,key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Hidden variables")
legend(x=0.9,y=1.1,c("TEA1","TEA2","TEA3"),title="Side Color",fill=tmp.pal[c("TEA1","TEA2","TEA3")],cex=0.6)
dev.off()

### Affinity matrix between hidden units
dm = as.matrix(distances(hid0[,hvar]))
am = exp(-dm^2/2)
pdf("Figures/hvar_affinity.pdf")
par(xpd=TRUE)
heatmap.2(am, trace="none", cexCol=0.5, col=colorRampPalette(c("white","red"))(32), scale="none",dendro="row",srtCol=270,offsetCol=0,offsetRow=0,cexRow=0.5,margin=c(3,3),RowSideColors=c1,ColSideColors=c1,density.info="none",key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Affinity")
legend(x=0.9,y=1.1,c("TEA1","TEA2","TEA3"),title="Side Color",fill=tmp.pal[c("TEA1","TEA2","TEA3")],cex=0.6)
dev.off()

### Correlation with TEA clusters
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

# Filter significantly correlated hidden units
hsig = which(abs(spcor)>0.65) # id of hidden units with significant correlation with TEA

pdf("Figures/spcorr.pdf", width=10, height=6)
par(mfrow=c(2,3))
par(mar=c(3,1,3,1))

for(i in hsig){
  if(i==36){plot.new()}
  boxplot(hid0[,i]~tea, main=paste("H",i,sep=""), cex=0.5, ylim=c(-0.2,1),col=tmp.pal[1:3])
  text(x=2, y=-0.11, paste("Spearman correlation = ",round(spcor[i],3), sep=""),cex=1)
}
dev.off()

### Hidden layers correlated with TEA clusters (Hsig)
### Heatmap and hierarchical clustering of Hsig
pdf("Figures/hsig.pdf")
par(xpd=TRUE)
heatmap.2(hid0[,hsig], trace="none", cexCol=1, col=colorRampPalette(c("white","red"))(32), scale="none",dendro="row",srtCol=270,offsetCol=0,offsetRow=0,adjCol=c(0.5,1),cexRow=0.5,margin=c(3,3),RowSideColors=c1,density.info="none",labRow=NA,key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Hidden variables")
legend(x=0.9,y=1.1,c("TEA1","TEA2","TEA3"),title="Side Color",fill=tmp.pal[c("TEA1","TEA2","TEA3")],cex=0.6)
dev.off()

### Affinity matrix between Hsig
pdf("Figures/hsig_affinity.pdf")
dm <- as.matrix(distances(hid0[,hsig]))
am <- exp(-dm^2/2)
par(xpd=TRUE)
heatmap.2(am, trace="none", cexCol=0.5, col=colorRampPalette(c("white","red"))(32), scale="none",dendro="row",srtCol=270,offsetCol=0,offsetRow=0,cexRow=0.5,margin=c(3,3),RowSideColors=c1,ColSideColors=c1,density.info="none",key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Affinity")
legend(x=0.9,y=1.1,c("TEA1","TEA2","TEA3"),title="Side Color",fill=tmp.pal[c("TEA1","TEA2","TEA3")],cex=0.6)
dev.off()

### Weight distribution of Hsig
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


### Load annotation and gene names
myinf2 = "data/GPL6244.annot"

### Read in the gene name for each probeset
conIn = file(myinf2, "r")
data = readLines(conIn, -1)
close(conIn)
mysta = grep("platform_table_begin", data)
myend = grep("platform_table_end", data)
data = read.table(myinf2, sep="\t", header=T, row.names=1, skip=mysta, nrows=myend-mysta-2, quote="", comment.char ="")
tmp = row.names(data)

GPL = data[, "GenBank.Accession"]
GS = data[, "Gene.symbol"]
GT = data[,"Gene.title"]
names(GPL) = tmp
names(GS) = tmp
names(GT) = tmp

gene.id = as.character(GPL[probe.id])

### Map probe ids to gene names
probes.gene.names = unlist(lapply(as.character(GS[probe.id]), function(x){strsplit(x, "/")[[1]][1]}))
names(probes.gene.names) = probe.id

save.image("data.RData")
