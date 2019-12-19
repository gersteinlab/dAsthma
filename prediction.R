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
library(reshape2)

################## Correlation with clinical traits
### All
clinical <- read.table("data/clinical.csv",sep=",")
colnames(clinical) = c("Array_ID", "PRED.FEV1", "PRED.FVC", "POST.FEV1", "POST.FVC", "POST.FEV1/FVC", "PRE.FEV1", "PRE.FVC", "PRE.FEV1/FVC")
rownames(clinical) = clinical[,1]
clinical = clinical[,-1]

clinical.spcor = data.frame(matrix(0L, nrow=6, ncol=length(hvar)))
rownames(clinical.spcor) = colnames(clinical)[3:8]
colnames(clinical.spcor) = paste("H", hvar,sep="")

for(i in seq(1,length(hvar))){
  for(j in seq(1,6)){
    hv = hid0[,hvar[i]]
    t = as.numeric(as.character(clinical[non.ctrl,j+2]))
    id = which(!is.na(t))
    clinical.spcor[j,i] = cor(hv[id], t[id], method="spearman")
  }
}

pdf("Figures/clinical_spcorr.pdf",width=6, height=10)
heatmap.2(t(clinical.spcor), trace="none", col=colorRampPalette(c("white","steelblue"))(32), Rowv="none", Colv="none", dendrogram="none",colsep=1:nrow(clinical.spcor), rowsep=1:ncol(clinical.spcor),sepcolor="white")
dev.off()

### Ratio
pdf("Figures/ratio_spcorr.pdf",width=16, height=6)
barplot(as.matrix(clinical.spcor[c(6,3),])[,order(clinical.spcor[3,], decreasing=T)],beside=TRUE, col=rep(c("navy","red"),length(hsig)), cex.axis = 0.6)
legend("bottomleft", c("PRE.FEV1/FVC","POST.FEV1/FVC"), fill=c("navy","red"), col=c("navy","red"), cex=0.6, pt.cex = 1)
dev.off()

################## Regression analysis for clinical traits
variance_explained <- function(y_true, y_pred){
  y_diff_var = var(y_true - y_pred)
  y_true_var = var(y_true)
  return(1-y_diff_var/y_true_var)
}

### Linear regression
reg_performance <- function(xx, yy){
  y = yy[which(!is.na(yy))]
  x = xx[which(!is.na(yy)),]
  
  N = length(y)
  mse_all = rep(0,100)
  ve_all = rep(0,100)
  for(n in seq(1,100)){
    id.parts = sample(1:4,size=N,replace=TRUE)
    mse = rep(0,4)
    ve = rep(0,4)
    for(i in seq(1,4)){
      id = which(id.parts==i)
      dat.train = data.frame(x = x[-id,], y=y[-id])
      dat.test = data.frame(x = x[id,], y=y[id])
      
      lr = lm(y~., data=dat.train)
      y.pred = predict(lr, dat.test)
      mse[i] = mean((y.pred-y[id])**2)
      ve[i] = variance_explained(y[id], y.pred)
    }
    mse_all[n] = mean(mse)
    ve_all[n] = mean(ve)
    
  }
  print(paste("MSE =", mean(mse_all)))
  print(paste("Variance explained =", mean(ve_all)))
  return(mean(mse_all))
}

### SVM
svm_performance <- function(xx, yy){
  y = yy[which(!is.na(yy))]
  x = xx[which(!is.na(yy)),]
  
  N = length(y)
  mse_all = rep(0,100)
  ve_all = rep(0,100)
  for(n in seq(1,100)){
    id.parts = sample(1:4,size=N,replace=TRUE)
    mse = rep(0,4)
    ve = rep(0,4)
    for(i in seq(1,4)){
      id = which(id.parts==i)
      dat.train = data.frame(x = x[-id,], y=y[-id])
      dat.test = data.frame(x = x[id,], y=y[id])
      
      lr = svm(y~., data=dat.train)
      y.pred = predict(lr, dat.test)
      mse[i] = mean((y.pred-y[id])**2)
      ve[i] = variance_explained(y[id], y.pred)
    }
    mse_all[n] = mean(mse)
    ve_all[n] = mean(ve)
  }
  print(paste("MSE =", mean(mse_all)))
  print(paste("Variance explained =", mean(ve_all)))
  return(mean(mse_all))
}

### Lasso
lasso_performance <- function(xx, yy, lambda=10^-2){
  y = yy[which(!is.na(yy))]
  x = xx[which(!is.na(yy)),]
  
  N = length(y)
  mse_all = rep(0,100)
  ve_all = rep(0,100)
  for(n in seq(1,100)){
    id.parts = sample(1:4,size=N,replace=TRUE)
    mse = rep(0,4)
    ve = rep(0,4)
    for(i in seq(1,4)){
      id = which(id.parts==i)
      dat.train = data.frame(x = x[-id,], y=y[-id])
      dat.test = data.frame(x = x[id,], y=y[id])
      
      lr = glmnet(x[-id,], y[-id], alpha = 1, lambda = lambda)
      y.pred = predict(lr, x[id,])
      mse[i] = mean((y.pred-y[id])**2)
      ve[i] = variance_explained(y[id], y.pred)
    }
    mse_all[n] = mean(mse)
    ve_all[n] = mean(ve)
  }
  print(paste("MSE =", mean(mse_all)))
  print(paste("Variance explained =", mean(ve_all)))
  return(mean(mse_all))
}

post.ratio = as.numeric(as.character(clinical$`POST.FEV1/FVC`))[non.ctrl]
reg_performance(h.sig, post.ratio)
reg_performance(h.var, post.ratio)

svm_performance(h.sig, post.ratio)
svm_performance(h.var, post.ratio)
svm_performance(x, post.ratio)
svm_performance(x[,imp.probes], post.ratio)

lambda = 10^-2
lasso_performance(h.sig, post.ratio, lambda)
lasso_performance(h.var, post.ratio, lambda)
lasso_performance(x, post.ratio, lambda)
lasso_performance(x[,imp.probes], post.ratio, lambda)

pre.ratio = as.numeric(as.character(clinical$`PRE.FEV1/FVC`))[non.ctrl]
reg_performance(h.sig, pre.ratio)
reg_performance(h.var, pre.ratio)

svm_performance(h.sig, pre.ratio)
svm_performance(h.var, pre.ratio)
svm_performance(x, pre.ratio)
svm_performance(x[,imp.probes], pre.ratio)

lambda = 10^-2
lasso_performance(h.sig, pre.ratio, lambda)
lasso_performance(h.var, pre.ratio, lambda)
lasso_performance(x, pre.ratio, lambda)
lasso_performance(x[,imp.probes], pre.ratio, lambda)

### DEG
tea1_deg_probe <- as.character(read.table("data/tea_degs/tea1_degs_id.txt")$V1)
tea3_deg_probe <- as.character(read.table("data/tea_degs/tea3_degs_id.txt")$V1)

deg = unique(c(tea1_deg_probe, tea3_deg_probe))
svm_performance(t(dat[deg, ]), post.ratio)
lasso_performance(t(dat[deg, ]), post.ratio, lambda)
svm_performance(t(dat[deg, ]), pre.ratio)
lasso_performance(t(dat[deg, ]), pre.ratio, lambda)



### Plot
plot_svm_performance <- function(xx, yy, path=NULL, main=NULL){
  if(!is.null(path)){pdf(path)}
  y = yy[which(!is.na(yy))]
  x = xx[which(!is.na(yy)),]
  
  N = length(y)
  id.parts = sample(1:4,size=N,replace=TRUE)
  dat.all = data.frame(x=x,y=y)

  id = which(id.parts==1)
  dat.train = dat.all[-id,]
  dat.test = dat.all[id,]

  lr = svm(y~., data=dat.train)
  y.pred = predict(lr, dat.all)
  ymin = min(c(y.pred, y))
  ymax = max(c(y.pred, y))

  plot(x=y, y=y.pred, xlab="True", ylab="Predicted", xlim=c(ymin, ymax), ylim=c(ymin,ymax), col="blue", main=main)
  abline(lm(y.pred~y), col="red")
  lines(x=c(0,1), y=c(0,1), col="black", lty=3)
}

plot_lasso_performance <- function(xx, yy, path=NULL, main=NULL){
  if(!is.null(path)){pdf(path)}
  y = yy[which(!is.na(yy))]
  x = xx[which(!is.na(yy)),]
  
  N = length(y)
  id.parts = sample(1:4,size=N,replace=TRUE)
  dat.all = data.frame(x=x,y=y)
  
  id = which(id.parts==1)
  dat.train = dat.all[-id,]
  dat.test = dat.all[id,]
  
  lr = glmnet(x[-id,], y[-id], alpha = 1, lambda = lambda)
  y.pred = predict(lr, x)
  ymin = min(c(y.pred, y))
  ymax = max(c(y.pred, y))
  
  plot(x=y, y=y.pred, xlab="True", ylab="Predicted", xlim=c(ymin, ymax), ylim=c(ymin,ymax), col="blue", main=main)
  abline(lm(y.pred~y), col="red")
  lines(x=c(0,1), y=c(0,1), col="black", lty=3)
}

svm_corr <- function(xx, yy, path=NULL, main=NULL){
  if(!is.null(path)){pdf(path)}
  y = yy[which(!is.na(yy))]
  x = xx[which(!is.na(yy)),]
  
  N = length(y)
  id.parts = sample(1:4,size=N,replace=TRUE)
  dat.all = data.frame(x=x,y=y)
  
  y.pred.list = NULL
  y.true.list = NULL
  for(i in seq(1,4)){
    id = which(id.parts==1)
    dat.train = dat.all[-id,]
    dat.test = dat.all[id,]
    
    lr = svm(y~., data=dat.train)
    y.pred = predict(lr, dat.test)
    
    y.pred.list = c(y.pred.list, y.pred)
    y.true.list = c(y.true.list, dat.test$y)
    
  }
  
  return(cor(y.pred.list, y.true.list))
}

lasso_corr <- function(xx, yy, path=NULL, main=NULL){
  if(!is.null(path)){pdf(path)}
  y = yy[which(!is.na(yy))]
  x = xx[which(!is.na(yy)),]
  
  N = length(y)
  id.parts = sample(1:4,size=N,replace=TRUE)
  dat.all = data.frame(x=x,y=y)
  
  id = which(id.parts==1)
  dat.train = dat.all[-id,]
  dat.test = dat.all[id,]
  
  lr = glmnet(x[-id,], y[-id], alpha = 1, lambda = lambda)
  y.pred = predict(lr, x[id,])
  return(cor(y.pred, y[id]))
}

plot_svm_performance(x[,imp.probes], pre.ratio, path="Figures/pred.pre.ratio.svm.pdf", main="PRE.FEV1/FVC (SVM)")
dev.off()

plot_svm_performance(x[,imp.probes], post.ratio, path="Figures/pred.post.ratio.svm.pdf", main="POST.FEV1/FVC (SVM)")
dev.off()

plot_lasso_performance(x[,imp.probes], pre.ratio, path="Figures/pred.pre.ratio.lasso.pdf", main="PRE.FEV1/FVC (LASSO)")
dev.off()

plot_lasso_performance(x[,imp.probes], post.ratio, path="Figures/pred.post.ratio.lasso.pdf", main="POST.FEV1/FVC (LASSO)")
dev.off()


### Plot bar plot
SVR.mse <- read.table(textConnection("Traits  H	Hsig	top.weighted	top.weighted.selected
POST.FEV1/FVC	0.0151	0.0138	0.0119	0.0116
PRE.FEV1/FVC	0.0166	0.0147  0.0122	0.0122"), header=T, row.names=1)

LASSO.mse <- read.table(textConnection("Traits  H	Hsig	top.weighted	top.weighted.selected
POST.FEV1/FVC	0.0139	0.0135	0.0128	0.0121
PRE.FEV1/FVC  0.0141	0.0139	0.0138	0.0130"), header=T, row.names=1)

pre.mse = data.frame(mse = c(t(SVR.mse[2,]), t(LASSO.mse[2,])),
                     method = c(rep("SVR",4), rep("LASSO",4)),
                     data = c(colnames(SVR.mse), colnames(LASSO.mse)))

post.mse = data.frame(mse = -c(t(SVR.mse[1,]), t(LASSO.mse[1,])),
                      method = c(rep("SVR",4), rep("LASSO",4)),
                      data = c(colnames(SVR.mse), colnames(LASSO.mse)))

ggplot(mapping=aes(x=data, y=mse, fill=method, group=method)) + 
  geom_bar(data=pre.mse, stat="identity", color="black", position=position_dodge())+
  geom_bar(data=post.mse, stat="identity", color="black", position=position_dodge())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
