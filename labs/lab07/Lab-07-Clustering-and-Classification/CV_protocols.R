library(kernlab)
library(ROCR)
library(randomForest)

svm.kfoldAUC = function(x, y, k=5, ...){
  n=length(y)
  ### Attempt a balanced split; only possible if k is smaller 
  ### than the number of cases in the smallest class
  split.idx = split(unlist(tapply(1:n, y, sample)), rep(1:k, length=n))
  aucs = rep(NA, k)
  split.sample = list()
  x[!is.finite(x)] = 0
  y = factor(y)
  for(i in 1:k){
    sps = table(1:n %in% split.idx[[i]], y)
    if(!any(sps==0)){
      svp <- ksvm(x[-split.idx[[i]],], 
                  y[-split.idx[[i]]], 
                  type="C-svc", 
                  kernel='vanilladot', 
                  C=1, scaled=F, kpar=list(), ...)
      aucs[i] = 
        performance(
          prediction(
            predict(svp, x[split.idx[[i]],], type='decision'), 
            y[split.idx[[i]]]), 
          measure='auc')@y.values[[1]]
    }
    
    split.sample = append(split.sample, sps)
  }
  list(aucs=aucs, splits = split.idx, splits.sample.size = split.sample)
}


rf.kfoldAUC = function(x, y, k=5, ...){
  n=length(y)
  split.idx = split(unlist(tapply(1:n, y, sample)), rep(1:k, length=n))
  aucs = rep(NA, k)
  split.sample = list()
  x[!is.finite(x)] = 0
  y = factor(y)
  
  for(i in 1:k){
    sps = table(1:n %in% split.idx[[i]], y)
    if(!any(sps==0)){
      rfm <- randomForest(x[-split.idx[[i]],], 
                          y[-split.idx[[i]]], ...)
      aucs[i] = performance(
        prediction(predict(rfm, x[split.idx[[i]],], type="prob")[,2],
                   y[split.idx[[i]]]), 
        measure='auc')@y.values[[1]]
    }
    
    split.sample = append(split.sample, sps)
  }
  list(aucs=aucs, splits = split.idx, splits.sample.size = split.sample)  
}