require(gplots)
require(RColorBrewer)
require(distances)
require(preprocessCore)

dat.raw <- read.table("data/sputum_combat_RIN_adjusted_all.txt",sep="\t")
dat <- normalize.quantiles(as.matrix(dat.raw))
rownames(dat) <- rownames(dat.raw)
colnames(dat) <- colnames(dat.raw)

frq <- read.table("data/freq.sample.clinic.txt")
non.ctrl=which(frq$V9!="CONTROL")
sev <- as.factor(frq$V9[non.ctrl])
sev_sub_id <- c(which(sev=="MILD"), which(sev=="SEVERE"))
sev_sub <- sev[sev_sub_id]

hid <- as.matrix(read.table("data/New.sput.cycle100.h50.lr0.1.err0.01.hidden.txt",sep="\t"))
wt <- as.matrix(read.table("data/New.sput.cycle100.h50.lr0.1.err0.01.weight.txt",sep="\t"))
wt <- wt[,-1]
probe.id <- rownames(dat.raw)
colnames(wt) <- probe.id

tea <- read.table("data/TEA_kegg.dist_kmeans_3.txt")
tea <- as.factor(paste("tea",tea$V2, sep=""))

hid0 = hid[non.ctrl,-1]
colnames(hid0) = paste("V",seq(1,50),sep="")
hvar = which(apply(hid0, 2, var)>0.001) # id of hidden units with variance
hid0 = hid0[sev_sub_id,] # SEVERE and MILD samples
  
tmp.pal=brewer.pal(5,"Set2");
names(tmp.pal)<-c("CONTROL","MILD","MODERATE","SEVERE");

c1 <- tmp.pal[as.character(sev_sub)];

### Heatmap
par(xpd=TRUE)
heatmap.2(hid0[,hvar], trace="none", cexCol=1, col=colorRampPalette(c("white","red"))(256), scale="none",dendro="row",srtCol=270,offsetCol=0,offsetRow=0,adjCol=c(0.5,1),cexRow=0.5,margin=c(3,3),RowSideColors=c1,density.info="none",labRow=NA,key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Hidden variables")
legend(x=0.9,y=1.1,c("MILD","SEVERE"),title="Side Color",fill=tmp.pal[c("MILD","SEVERE")],cex=0.6)

### Affinity matrix
dm <- as.matrix(distances(hid0[,hvar]))
am <- exp(-dm^2/2)
par(xpd=TRUE)
heatmap.2(am, trace="none", cexCol=0.5, col=colorRampPalette(c("white","red"))(256), scale="none",dendro="row",srtCol=270,offsetCol=0,offsetRow=0,cexRow=0.5,margin=c(3,3),RowSideColors=c1,ColSideColors=c1,density.info="none",key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Affinity")
legend(x=0.9,y=1.1,c("MILD","SEVERE"),title="Side Color",fill=tmp.pal[c("MILD","SEVERE")],cex=0.6)

### Correlation with cluster
n = length(hid0[1,])

par(mfrow=c(4,6))
par(mar=c(3,1,3,1))

wilcox.pval = rep(1, n)

for(i in hvar){
  hi = hid0[,i]
  boxplot(hid0[,i]~sev_sub, main=paste("V",i,sep=""), cex=0.5, ylim=c(-0.2,1))
  wilcox.pval[i] = wilcox.test(hi[sev_sub=="MILD"], hi[sev_sub=="SEVERE"], exact = F)$p.value
  text(x=2.5, y=-0.11, "Wilcox pval =",cex=0.8)
  text(x=2.5, y=-0.2, wilcox.pval[i], cex=0.8)
}

hsig = which(wilcox.pval<0.05) # id of hidden units with significant correlation with tea

### Heatmap
par(xpd=TRUE)
heatmap.2(hid0[,hsig], trace="none", cexCol=1, col=colorRampPalette(c("white","red"))(256), scale="none",dendro="row",srtCol=270,offsetCol=0,offsetRow=0,adjCol=c(0.5,1),cexRow=0.5,margin=c(3,3),RowSideColors=c1,density.info="none",labRow=NA,key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Hidden variables")
legend(x=0.9,y=1.1,c("MILD","SEVERE"),title="Side Color",fill=tmp.pal[c("MILD","SEVERE")],cex=0.6)

### Affinity matrix
dm <- as.matrix(distances(hid0[,hsig]))
am <- exp(-dm^2/2)
par(xpd=TRUE)
heatmap.2(am, trace="none", cexCol=0.5, col=colorRampPalette(c("white","red"))(256), scale="none",dendro="row",srtCol=270,offsetCol=0,offsetRow=0,cexRow=0.5,margin=c(3,3),RowSideColors=c1,ColSideColors=c1,density.info="none",key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Affinity")
legend(x=0.9,y=1.1,c("MILD","SEVERE"),title="Side Color",fill=tmp.pal[c("MILD","SEVERE")],cex=0.6)


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

GPL = data[, "Gene.symbol"]
#GPL = data[, "GenBank.Accession"]
names(GPL) = tmp

gene.id = as.character(GPL[probe.id])


### Top scored genes
top.gene = list()
top.probe = list()

top.gene.GO = list()
N = length(probe.id)


#for(i in hsig){
#  v.id = paste("V",i,sep="")
#  top.id = order(wt[i,])[(N-100+1):N]
#  tmpv = rep(0,N)
#  names(tmpv) = probe.id
#  tmpv[top.id] = 1
#  tmpv = as.factor(tmpv)
#  top.100[[v.id]] = tmpv
#  GOdata <- new("topGOdata",
#                description = "topGO", ontology = "BP",
#                allGenes = tmpv, geneSel = NULL,
#                nodeSize = 10,
#                annot = annFUN.db, affyLib = "hugene10sttranscriptcluster")
#  resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher") 
#  sig.tab <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 20, numChar=1000)

#  top.100.GO[[v.id]] = sig.tab[,1:2]
#}



for(i in hsig){
  v.id = paste("V",i,sep="")
  wv = wt[i,]
  meanv = mean(wv)
  sdv = sd(wv)
  top.id = which(abs(wv-meanv)/sdv>=3)
  top.gene[[v.id]] = gene.id[top.id]
  top.probe[[v.id]] = probe.id[top.id]
  #write(gene.id[top.id], file=paste("/Users/tianxiao/Documents/Lab/Asthma/V",i,".txt",sep=""))
}


### Heatmap of selected genes
probe.list = unique(unlist(top.probe))

par(xpd=TRUE)
heatmap.2(as.matrix(dat[probe.list,sev_sub_id]), trace="none", col=colorRampPalette(c("green","black","red"))(256), scale="row",offsetCol=0,offsetRow=0,adjCol=c(0.5,1),margin=c(3,3),ColSideColors=c1,density.info="none",labRow=NA, labCol=NA,key.title="",key.xlab="",key.ylab="",key.par=list(cex=0.5,cex.main=0.5,cex.lab=0.2), lhei = c(1, 4),lwid=c(1,4),keysize=1, main="Expression")
legend(x=0.9,y=1.1,c("MILD","SEVERE"),title="Side Color",fill=tmp.pal[c("MILD","SEVERE")],cex=0.6)


