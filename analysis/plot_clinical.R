### Correlation between clinical information
### Plot correlation with clinical traits
load('data.RData')

clinical_all = read.table("data/clinical_all.csv",sep=",", header=T)
rownames(clinical_all) = clinical_all[,1]
clinical_all = clinical_all[,-1]

na.count = apply(clinical_all, 2, function(x){sum(is.na(x))})
clinical_filtered = clinical_all[,na.count<=10]
clinical_filtered = clinical_filtered[,-c(43,57,58,59,60,62)]


### Correlation between H and clinical traits
clinical_all.spcor = data.frame(matrix(0L, nrow=ncol(clinical_filtered), ncol=length(hvar)))
rownames(clinical_all.spcor) = colnames(clinical_filtered)
colnames(clinical_all.spcor) = paste("H", hvar,sep="")

for(i in seq(1,length(hvar))){
  for(j in seq(1,ncol(clinical_filtered))){
    hv = hid0[,hvar[i]]
    t = clinical_filtered[non.ctrl,j]
    id = which(!is.na(t))
    t = as.numeric(t)
    clinical_all.spcor[j,i] = cor(hv[id], t[id], method="spearman")
  }
}

pdf("Figures/H_clinical_spcorr.pdf",width=6, height=12)
heatmap.2(t(t(clinical_all.spcor)), trace="none", col=colorRampPalette(c("blue","white","red"))(32), 
          Rowv="none", dendrogram="none",
          colsep=1:ncol(clinical_all.spcor), rowsep=1:nrow(clinical_all.spcor),sepcolor="white")
dev.off()

### Correlation between clinical traits
clinical_all.inter.spcor = data.frame(matrix(0L, ncol=ncol(clinical_filtered), nrow=ncol(clinical_filtered)))
rownames(clinical_all.inter.spcor) = colnames(clinical_filtered)
colnames(clinical_all.inter.spcor) = colnames(clinical_filtered)
for(i in seq(1,ncol(clinical_filtered))){
  for(j in seq(1,ncol(clinical_filtered))){
    ti = clinical_filtered[non.ctrl,i];tj = clinical_filtered[non.ctrl,j]
    id = which((!is.na(ti))*(!is.na(tj))>0)
    ti = as.numeric(ti);tj = as.numeric(tj)
    clinical_all.inter.spcor[j,i] = cor(ti[id], tj[id], method="spearman")
  }
}

pdf("Figures/clinical_inter_spcorr.pdf",width=12, height=12)
heatmap.2(t(t(clinical_all.inter.spcor)), trace="none", col=colorRampPalette(c("blue","white","red"))(32),
          dendrogram="none",colsep=1:ncol(clinical_all.inter.spcor), rowsep=1:nrow(clinical_all.inter.spcor),sepcolor="white")
dev.off()

### Boxplot and scatter plot
par(mfrow=c(2,3))

for(i in hsig){
  if(i==36){plot.new()}
  boxplot(as.numeric(hid0[,i])~as.character(clinical_filtered$Gender[non.ctrl]), main=paste("H",i, sep=""))
}


par(mfrow=c(2,3))

for(i in hsig){
  if(i==36){plot.new()}
  plot(as.numeric(hid0[,i])~as.character(clinical_filtered$BDR[non.ctrl]), main=paste("H",i, sep=""))
}

