### Feature selection and prediction
### Prediction of TEA clusters and asthma severity based on expression of selected genes/hidden unit values

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

###### Filter genes for prediction
### Top weighted genes for each Hsig
probes.ribosomal = read.table("data/probes_ribosomal.txt")$V1
idx.non.ribosomal = (probe.id %in% probes.ribosomal == FALSE)

wt.non.ribosomal = wt[,idx.non.ribosomal]

# Top and bottom-weighted genes in Hsig
top.gene = list()
top.probe = list()
top.gene.GO = list()

bot.gene = list()
bot.probe = list()
bot.gene.GO = list()

N = length(probe.id)
for(i in hsig){
  v.id = paste("V",i,sep="")
  wv = wt.non.ribosomal[i,]
  names(wv) = probes.gene.names[idx.non.ribosomal]
  #wv = wt[i,var.probe.idx]
  meanv = mean(wv)
  sdv = sd(wv)
  
  top.id = order(wv, decreasing = T)[1:200]
  top.probe[[v.id]] = probe.id[idx.non.ribosomal][top.id]
  top.gene[[v.id]] = probes.gene.names[idx.non.ribosomal][top.id]
  
  bot.id = order(wv)[1:200]
  bot.probe[[v.id]] = probe.id[idx.non.ribosomal][bot.id]
  bot.gene[[v.id]] = probes.gene.names[idx.non.ribosomal][bot.id]
  
  top.gene[[v.id]] = gene.id[top.id]
  write(gene.id[top.id], file=paste("genes_hidden_units/H",i,"_top.txt",sep=""))
  write(gene.id[bot.id], file=paste("genes_hidden_units/H",i,"_bot.txt",sep=""))
}

top.probe.list = unique(unlist(top.probe))
bot.probe.list = unique(unlist(bot.probe))
probe.list = unique(c(top.probe.list,bot.probe.list))

# Prepare data
x = t(dat[top.probe.list,non.ctrl]) # expression of top-weighted genes
h.sig = hid0[, hsig] # Hsig values
h.var = hid0[, hvar] # Hvar values

###################### Random forest prediction ######################
# Train the random forest classifier
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

# Plot random forest performance (ROC)
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

###################### Prediction of TEA clusters ######################
### TEA-related DEGs
tea1_deg_probe <- as.character(read.table("data/tea_degs/tea1_degs_id.txt")$V1)
tea3_deg_probe <- as.character(read.table("data/tea_degs/tea3_degs_id.txt")$V1)
tea_deg = unique(c(tea1_deg_probe, tea3_deg_probe))

tea_1and3_id = c(which(tea=="TEA1"), which(tea=="TEA3"))
t.tea = as.factor(as.numeric(tea[tea_1and3_id])==3)

rf_performance(h.sig[tea_1and3_id, ], t.tea)
rf_performance(h.var[tea_1and3_id, ], t.tea)
rf_performance(x[tea_1and3_id, ], t.tea)

###################### Prediction of severity ######################
sev_sub_id = c(which(sev=="MILD"), which(sev=="SEVERE"))
sev_sub = as.factor(as.character(sev[sev_sub_id]))

t.sev =  as.factor(as.numeric(sev[sev_sub_id])==3)

### Select severity-related genes based on random forest importance
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


N.feat = length(x[1,])
mean.imp = mean(gene.imp.mean)
std.imp = sd(gene.imp.mean)
#imp.probes = names(gene.imp.mean)[(gene.imp.mean-mean.imp)/std.imp > 1]
imp.probes = rownames(gene.imp)[order(gene.imp.mean)[(N.feat-49):N.feat]]
top.imp.probes = rownames(gene.imp)[order(gene.imp.mean)[(N.feat-9):N.feat]]

# Plot importance
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

### Select severity-related DEGs
group_deg_top <- function(dat, label, c1, c2, cutoff=50){
  names = rownames(dat)
  pval = c()
  dat1 = dat[,(label==c1)]
  dat2 = dat[,(label==c2)]
  for(i in seq(1, nrow(dat))){
    exp1 = dat1[i,]
    exp2 = dat2[i,]
    tt = t.test(exp1, exp2)
    pval = c(pval, tt$p.value)
  }
  pval = order(pval, decreasing = FALSE)
  deg = names[pval[1:cutoff]]
  return(deg)
}

sev_deg = group_deg_top(dat, sev, 'MILD', 'SEVERE', cutoff=50)

### Compare performance
rf_performance(t(dat[sev_deg, sev_sub_id]), t.sev)
rf_performance(x[sev_sub_id, imp.probes], t.sev)
rf_performance(h.sig[sev_sub_id, ], t.sev)
rf_performance(h.var[sev_sub_id, ], t.sev)
rf_performance(x[sev_sub_id, ], t.sev)

# Compare with random
rf_rand = rep(0,10)
for(i in seq(1,10)){
  rand.probes = sample(colnames(x),50)
  rf_rand[i] = rf_performance(x[sev_sub_id, rand.probes], t.sev)
}
mean(rf_rand)


### Plot ROC
plot_rf_performance(t(dat[tea_deg, sev_sub_id]), t.sev, path="roc.pdf", col="blue")
plot_rf_performance(t(dat[sev_deg, sev_sub_id]), t.sev, col="yellow")
plot_rf_performance(h.sig[sev_sub_id, ], t.sev, col="purple", add=T)
plot_rf_performance(h.var[sev_sub_id, ], t.sev, col="green", add=T)
plot_rf_performance(x[sev_sub_id, ], t.sev, col="cyan", add=T)
plot_rf_performance(x[sev_sub_id, imp.probes], t.sev, col="red", add=T)
dev.off()

save.image("data_pred.RData")
