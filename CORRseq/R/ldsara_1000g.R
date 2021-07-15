
#using excavtor2 normalization procedure and CODEX2 QC creteria 
for(chr in 1:22){
 setwd("/Users/xizhiluo/Dropbox/ldsara_seq/")
 # setwd("C:\\Users\\Adam Luo\\Dropbox\\ldsara_seq\\")

  library(modSaRa)   
  source("CNVcluster.R")  
  source("CNVout.R")  
  source("gausianMixture.R")  
  source("modifiedSaRa.R")
  source("multiSaRa.R")
  source("getOneBIC.R") 
  map=get(load(file=paste("/Volumes/Samsung USB/codex_ldsara/mapchr",chr,".RData",sep="")))
  lrr=get(load(paste("/Volumes/Samsung USB/codex_ldsara/lrrchr",chr,".RData",sep="")))
  
  #map=get(load(file=paste("C:\\Users\\Adam Luo\\Desktop\\codex_ldsara\\mapchr",chr,".RData",sep="")))
  #lrr=get(load(paste("C:\\Users\\Adam Luo\\Desktop\\codex_ldsara\\lrrchr",chr,".RData",sep="")))

  h=5
  h_vec=c(rep(1/h,h),rep(-1/h,h))
  empvar=vector()
  for(d in (h+1):(dim(lrr)[1]-h)){
    local_cov=cov(t(lrr[(d-h):(d+h-1),]))
    ecnvar=t(h_vec)%*%local_cov%*%h_vec
    empvar[d]=ecnvar
  }
  empsd1=sqrt(empvar[!is.na(empvar)])
  
  
  h=10
  h_vec=c(rep(1/h,h),rep(-1/h,h))
  empvar=vector()
  for(d in (h+1):(dim(lrr)[1]-h)){
    local_cov=cov(t(lrr[(d-h):(d+h-1),]))
    ecnvar=t(h_vec)%*%local_cov%*%h_vec
    empvar[d]=ecnvar
  }
  empsd2=sqrt(empvar[!is.na(empvar)])
  
  
  h=15
  h_vec=c(rep(1/h,h),rep(-1/h,h))
  empvar=vector()
  for(d in (h+1):(dim(lrr)[1]-h)){
    local_cov=cov(t(lrr[(d-h):(d+h-1),]))
    ecnvar=t(h_vec)%*%local_cov%*%h_vec
    empvar[d]=ecnvar
  }
  empsd3=sqrt(empvar[!is.na(empvar)])
  
  
  
  
  fInverse <-
    function(n = 10000, h= 10, hh = 2 * h, precise = 10000, emp,simT = 2000){        #need to be faster
      empirical = NULL
      for (i in 1 : simT){
        Y=vector()
        for(p in 1:n){
          Y[p]=rnorm(1,mean = 0,sd=emp[p]) 
        }
        LDF  =  Y
        LDF.pos      =  LDF
        LDF.neg      =  -LDF
        index.pos = localMax(y = LDF.pos, span = hh)
        pV.pos   = 1 - 2 * pnorm(LDF.pos[index.pos] / (emp[index.pos]))
        index.neg = localMax(y = LDF.neg, span = hh)
        pV.neg   = 1 - 2 * pnorm(LDF.neg[index.neg] / (emp[index.neg]))
        index <- c(index.pos, index.neg)
        pv <- c(pV.pos, pV.neg)
        pv <- pv[order(index)]
        index <- sort(index)
        len <- length(index)
        rm.list <- NULL
        for (j in 1 : (len - 1)) {
          if(index[j] >= index[j + 1] - h){
            rm.list <- c(rm.list, j)
          }
        }
        if (length(rm.list) > 0) {
          pv <- pv[-rm.list]
        }
        empirical <- c(empirical, pv)
        if (length(empirical) > 10 * precise) break
      }
      return(quantile(empirical, probs = c(0 : precise) / precise))
    }
  
  FINV=list()
  FINV[[1]]=fInverse(n=length(empsd1),h=5,precise =10000,emp = empsd1)
  FINV[[2]]=fInverse(n=length(empsd2),h=10,precise =10000,emp = empsd2)
  FINV[[3]]=fInverse(n=length(empsd3),h=15,precise =10000,emp = empsd3)
  
  
  #FINV1=list()
  #FINV1[[1]]=fInverse(n=length(empsd1),h=5,precise =length(empsd1),emp = empsd1)
  #FINV1[[2]]=fInverse(n=length(empsd2),h=10,precise =length(empsd2),emp = empsd2)
  #FINV1[[3]]=fInverse(n=length(empsd3),h=15,precise =length(empsd3),emp = empsd3)
  
  
  for(alpha in c(0.001,0.01,0.05,0.1)){
    output=NULL
    #cp <-  CNVout(lrr=lrr,map=map, FINV=FINV, alpha = alpha,  outname = paste("C:\\Users\\Adam Luo\\Desktop\\ldsaraout\\ldoutchr",chr,alpha,".csv",sep=""))$cp
    cp <-  CNVout(lrr=lrr,map=map, FINV=FINV, alpha = alpha,  outname = paste("/Users/xizhiluo/Desktop/yhat_out//ldoutchr",chr,alpha,".csv",sep=""))$cp
  }
}


Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}
Modes(lrr[,1])