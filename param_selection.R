### Parameter selection (mtry)
train.acc = NULL
test.acc = NULL

for(n in seq(2,50,5)){
  xx=x[sev_sub_id, rownames(gene.imp)[order(gene.imp.mean)[(N.feat-n):N.feat]]]
  id.neg = which(as.numeric(yy)==1)
  id.pos = which(as.numeric(yy)==2)
  N.neg = sum(as.numeric(yy)==1)
  N.pos = sum(as.numeric(yy)==2)
  acc_all = rep(0,10)
  dat.all = data.frame(label=yy,xx)
  
  predicted.probs = NULL
  true.lables = NULL
  predicted.train.probs = NULL
  true.train.lables = NULL
  
  id.parts.neg = sample(1:4,size=N.neg,replace=TRUE)
  id.parts.pos = sample(1:4,size=N.pos,replace=TRUE)
  acc = rep(0,4)
  for(i in seq(1,4)){
    id = c(id.neg[which(id.parts.neg==i)],id.pos[which(id.parts.pos==i)])
    dat.train = dat.all[-id,]
    dat.test = dat.all[id,]
    
    tunegrid <- expand.grid(.mtry=mtry)
    rf = train(label~., dat.train, tuneGrid=tunegrid)
    
    predicted.prob = predict(rf, dat.test, type="prob")[,"TRUE"]
    predicted.probs = c(predicted.probs, predicted.prob)
    true.lables = c(true.lables, dat.test$label)
    
    predicted.train.prob = predict(rf, dat.train, type="prob")[,"TRUE"]
    predicted.train.probs = c(predicted.train.probs, predicted.train.prob)
    true.train.lables = c(true.train.lables, dat.train$label)
    
    #result.roc = roc(dat.test$label, predicted.prob[,"TRUE"]) # Draw ROC curve.
    #cm = confusionMatrix(y[id], predict(rf, x[id,]))
    #acc[i] = cm$overall[1]
    #acc[i] = pROC::auc(result.roc)
  }
  result.roc = roc(true.lables, predicted.probs) # Draw ROC curve.
  train.roc = roc(true.train.lables, predicted.train.probs) # Draw ROC curve.
  acc_all = pROC::auc(result.roc)
  acc.train_all = pROC::auc(train.roc)
  
  train.acc = c(train.acc, acc.train_all)
  test.acc = c(test.acc, acc_all)
}


### DEG
xx = t(dat[,non.ctrl])

deg_pval = data.frame(gene=NULL, pval=NULL)

for(i in seq(1, ncol(xx))){
  dat.deg = data.frame(exp = xx[sev_sub_id,i], sev=t.sev)
  
  deg.t.test = t.test(exp~sev, data=dat.deg)

  deg_pval = rbind(deg_pval, data.frame(gene=colnames(xx)[i], 
                                        pval=deg.t.test$p.value))
}

deg.sort = deg_pval[order(deg_pval$pval),]
