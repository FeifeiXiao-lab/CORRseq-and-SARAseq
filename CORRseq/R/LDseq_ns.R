library(CODEX2)
library(modSaRa)
library(ExomeDepth)
library(lsa)
setwd("C:\\Users\\fqin\\Dropbox\\WES CNV\\calling results 3 methods\\Code\\ldsara_seq")
source("CNVcluster_new.R")
source("CNVout.R")  
source("gausianMixture.R")  
source("modifiedSaRa.R")
source("multiSaRa.R")
source("getOneBIC.R")
source("Normalization.R")
Path2ExomeRC <- "C:\\Users\\fqin\\Dropbox\\WES CNV\\EXCAVATOR PACKAGES\\lib\\R\\LibraryExomeRC.R"
source(Path2ExomeRC)
source("C:\\Users\\fqin\\Dropbox\\WES CNV\\calling results 3 methods\\LDseq\\Code\\exomeDepth select reference sets.R")

length="short"
cnv="del"
nstype="all"

fun_LD <- function(length ,cnv) {
    test_qc=get(load(file=paste("C:\\Users\\fqin\\Dropbox\\WES CNV\\calling results 3 methods\\CODEX2\\simulation\\RC.70S.",length,".",cnv,".chr1_2.RData",sep="")))
    wholedata=get(load(file=paste("C:\\Users\\fqin\\Dropbox\\WES CNV\\calling results 3 methods\\CODEX2\\simulation\\RC.93S.rmCNVs.chr1_2.RData",sep="")))
    sample_common <- get(load(file="C:\\Users\\fqin\\Dropbox\\WES CNV\\calling results 3 methods\\sample_common.RData"))
    ns <- wholedata[,-which(colnames(wholedata) %in% sample_common)]
    
    Y_qc <- cbind(test_qc,ns)
    ns_index <- colnames(wholedata)[-which(colnames(wholedata) %in% sample_common)]
    
    ref_qc=get(load(file=paste("C:\\Users\\fqin\\Dropbox\\WES CNV\\calling results 3 methods\\CODEX2\\simulation\\ref_dataframe.RData",sep="")))
    
    start=ref_qc$start
    end=ref_qc$end
    exon_size=end-start
    gc=ref_qc$gc
    map=ref_qc$mapp
    
    RCTL <- Y_qc/exon_size
    
    RC_norm=matrix(data=NA,nrow = dim(RCTL)[1],ncol = dim(RCTL)[2])
    rownames(RC_norm)=rownames(RCTL)
    colnames(RC_norm)=colnames(RCTL)
    for(sub in 1:dim(RCTL)[2]){
      ### Exon size normalization only if InTarget###
      step <- 5
      RCLNormListIn <- CorrectSize(RCTL[,sub],L=exon_size,step)
      RCLNormIn <- RCLNormListIn$RCNorm
      
      ### Mappability normalization ###
      step <- 0.01
      RCMAPNormListIn <- CorrectMAP(RCLNormIn,MAPContent=map,step)
      RCMAPNormIn <- RCMAPNormListIn$RCNorm
      
      
      ### GC-content Normalization ###
      step <- 5
      RCGCNormListIn <- CorrectGC(RCMAPNormIn,GCContent=gc,step)
      RCGCNormIn <- RCGCNormListIn$RCNorm
      
      RC_norm[,sub]=RCGCNormIn
    }


  for (nstype in c("all","random","correlation","cosine similarity","euclidean")){
    test_norm <- RC_norm[,!(colnames(RC_norm) %in% ns_index)]
    lrrdata=matrix(data=NA,nrow = dim(test_norm)[1],ncol = dim(test_norm)[2])
    rownames(lrrdata)=rownames(test_norm)
    colnames(lrrdata)=colnames(test_norm)
    
    if (nstype=="all"){ 
      message('Use all samples as reference set')
      ns_norm <- RC_norm[,ns_index]    
      for(i in 1:ncol(test_qc)){
        meanrow=apply(ns_norm,1,mean)
        meanrow[which(meanrow==0)] <- 1
        rc=RC_norm[,i]
        rc[which(rc==0)]=meanrow[which(rc==0)]
        lrrdata[,i]=log2(rc/meanrow)
      }
    }
    
    if (nstype=="random"){
      message('Random choosing samples as reference set')
      for(i in 1:ncol(test_qc)){
        set.seed(i)
        random <- sample(ns_index,round(length(ns_index)*0.4))
        ns_norm <- RC_norm[,ns_index] 
        meanrow=apply(ns_norm,1,mean)
        meanrow[which(meanrow==0)] <- 1
        rc=RC_norm[,i]
        rc[which(rc==0)]=meanrow[which(rc==0)]
        lrrdata[,i]=log2(rc/meanrow)
      }
    }

    if (nstype=="correlation"){
      message('Choose optimal refs using correction methods')
      ns_norm <- RC_norm[,ns_index]    
      for(i in 1:ncol(test_qc)){ 
        rc=RC_norm[,i]
        my.choice <- select.reference.set.new(method=nstype,
                                              test.counts = rc,
                                              sampname=colnames(RC_norm)[i],
                                              reference.counts = ns_norm,
                                              bin.length = exon_size/1000,
                                              n.bins.reduced = 3000)
        
        my.matrix <- as.matrix(ns_norm[, my.choice, drop = FALSE])
        my.reference.selected <- apply(X = my.matrix,
                                       MAR = 1,
                                       FUN = mean)
        my.reference.selected[which(my.reference.selected==0)] <- 1
        rc[which(rc==0)]=my.reference.selected[which(rc==0)]
        lrrdata[,i]=log2(rc/my.reference.selected)
      }
    }
    
    if (nstype=="cosine similarity"){
      message('Choose optimal refs using cosine similarity methods')
      ns_norm <- RC_norm[,ns_index]    
      for(i in 1:ncol(test_qc)){ 
        rc=RC_norm[,i]
        my.choice <- select.reference.set.new(method=nstype,
                                              test.counts = rc,
                                              sampname=colnames(RC_norm)[i],
                                              reference.counts = ns_norm,
                                              bin.length = exon_size/1000,
                                              n.bins.reduced = 3000)
        
        my.matrix <- as.matrix(ns_norm[, my.choice, drop = FALSE])
        my.reference.selected <- apply(X = my.matrix,
                                       MAR = 1,
                                       FUN = mean)
        my.reference.selected[which(my.reference.selected==0)] <- 1
        rc[which(rc==0)]=my.reference.selected[which(rc==0)]
        lrrdata[,i]=log2(rc/my.reference.selected)
      }
    }
    
    if (nstype=="euclidean"){
      message('Choose optimal refs using cosine euclidean distance')
      ns_norm <- RC_norm[,ns_index]    
      for(i in 1:ncol(test_qc)){ 
        rc=RC_norm[,i]
        my.choice <- select.reference.set.new(method=nstype,
                                              test.counts = rc,
                                              sampname=colnames(RC_norm)[i],
                                              reference.counts = ns_norm,
                                              bin.length = exon_size/1000,
                                              n.bins.reduced = 3000)
        
        my.matrix <- as.matrix(ns_norm[, my.choice, drop = FALSE])
        my.reference.selected <- apply(X = my.matrix,
                                       MAR = 1,
                                       FUN = mean)
        my.reference.selected[which(my.reference.selected==0)] <- 1
        rc[which(rc==0)]=my.reference.selected[which(rc==0)]
        lrrdata[,i]=log2(rc/my.reference.selected)
      }
    }
    
    
    log2R=lrrdata
    for(chr in 1:2){
      map=ref_qc[which(ref_qc$seqnames==paste("chr",chr,sep="")),1:3]
      lrr=log2R[which(ref_qc$seqnames==paste("chr",chr,sep="")),]
      
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
      
      
    
      fInverse <- function(n = 10000, h= 10, hh = 2 * h, precise = 10000, emp, simT = 2000){  #need to be faster
        empirical = NULL
        for (sim in 1 : simT){
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
      
      
      for (alpha in c(0.05)){
        output=NULL
        cp <-  CNVout(lrr=lrr,map=map, FINV=FINV, alpha = alpha, chr=chr,   
                      outname = paste("C:\\Users\\fqin\\Dropbox\\WES CNV\\calling results 3 methods\\LDseq\\Choose reference\\",
                                     cnv,".",length,"chr_",chr,"alpha_",alpha,"_",nstype,sep=""))$cp
      }
   }}
}



length.list <- c("short","medium","long")
cnv.list <- c("del","dup")

for (l in length.list) {
  for (c in cnv.list) {
        fun_LD(length=l ,cnv=c)
  }
}



fun_LD(length="short" ,cnv="dup", nstype="all")







### TPR FPR ####
### del dup are indentified in getting freq ###
setwd("C:\\Users\\fqin\\Dropbox\\WES CNV\\calling results 3 methods\\simu93s")

cpmedium <- get(load(file="cplong.RData"))
cpmedium <- get(load(file="cpmedium.RData"))
cpshort <- get(load("cpshort.RData"))



setwd("C:\\Users\\fqin\\Dropbox\\WES CNV\\calling results 3 methods")

bed.out <- get(load(file="bedout.RData"))
exon <- get(load(file="exons.new2.RData"))
exon$Chr <- paste0("chr",exon$chr)

sample_common <- get(load(file="C:\\Users\\fqin\\Dropbox\\WES CNV\\calling results 3 methods\\sample_common.RData"))
sample_common <- sample_common[order(sample_common)]
sample_common <- data.frame(sample=sample_common)
sample_common$sampleID <- 1:70


fun_TPFP <- function(name, CNVtype, cp, alpha, nstype){
  
  obsCNV.whole1 <- read.csv(paste0("C:\\Users\\fqin\\Dropbox\\WES CNV\\calling results 3 methods\\LDseq\\Choose reference\\",CNVtype,".",name,"chr_1alpha_",alpha,"_",nstype,".csv"),header=F)
  obsCNV.whole2 <- read.csv(paste0("C:\\Users\\fqin\\Dropbox\\WES CNV\\calling results 3 methods\\LDseq\\Choose reference\\",CNVtype,".",name,"chr_2alpha_",alpha,"_",nstype,".csv"),header=F)
  
  obsCNV.whole <- rbind(obsCNV.whole1,obsCNV.whole2)
  colnames(obsCNV.whole) <- c("sampleID","chr","start","end","width","nSNP","type")
  
  if (CNVtype=="del"){
    obsCNV <- obsCNV.whole[which(obsCNV.whole$type=="del"),]
  }
  if (CNVtype=="dup"){
    obsCNV <- obsCNV.whole[which(obsCNV.whole$type=="dup"),]
  }
  
  TrueCNV <- data.frame(chr=exon$chr[cp[1:60,1]],start=exon$start[cp[1:60,1]],end=exon$end[cp[1:60,2]])
  TrueCNVR=makeGRangesFromDataFrame(TrueCNV)
  
  length=0
  for (s in 1:70){
    obsCNV1 <- obsCNV[which(obsCNV$sample==s),]
    if (nrow(obsCNV1) > 0) {
      obsCNV1 <- obsCNV1[,c(2,3,4)]
      colnames(obsCNV1) <- c("chr","start","end")
      obsCNVR1 <- makeGRangesFromDataFrame(obsCNV1)
      
      ovCNV1 <- findOverlaps(obsCNVR1,TrueCNVR)
      cnvid1= unique(subjectHits(ovCNV1))
      length <- length+length(cnvid1)
    }
  }
  
  tp <- round(length/4200,2)
  
  obsCNV <- obsCNV[,c(2,3,4)]
  colnames(obsCNV) <- c("chr","start","end")
  obsCNVR <- makeGRangesFromDataFrame(obsCNV)
  
  ovCNV <- findOverlaps(obsCNVR,TrueCNVR)
  cnvid= unique(queryHits(ovCNV))
  
  fp <- (nrow(obsCNV.whole)-length(cnvid))/(nrow(exon)-4200)
  
  res <- data.frame(tp=tp,fp=fp)
  return(res)
}



length.list <- c("short")
cnv.list <- c("del","dup")
nstype.list <- c("all","random","correlation","cosine similarity","euclidean")
alphalist <- c(0.05,0.1,0.001,0.0001)



length.list <- c("short")
cnv.list <- c("dup")
nstype.list <- c("all")
alphalist <- c(0.05)


TPFP <- matrix(NA,1,6)
colnames(TPFP) <- c("Length","CNV","ref.Method","alpha","TP.rate","FP.rate")
i=0
for (l in length.list) {
  for (c in cnv.list) {
    for (m in nstype.list) {
      for (a in alphalist){
        if (l=="short"){
          res <- fun_TPFP(name=l,CNVtype=c,cp=cpshort, alpha=a, nstype=m)
        }
        if (l=="medium"){
          res <- fun_TPFP(name=l,CNVtype=c,cp=cpmedium, alpha=a, nstype=m)
        }    
        if (l=="long"){
          res <- fun_TPFP(name=l,CNVtype=c,cp=cplong, alpha=a, nstype=m)
        }
        i=i+1
        TPFP[i,1] <- l
        TPFP[i,2] <- c
        TPFP[i,3] <- m
        TPFP[i,4] <- a
        TPFP[i,5] <- res[,1]
        TPFP[i,6] <- res[,2]
      }
    }
  }
}

prop.table <- TPFP
prop.table



setwd("C:\\Users\\fqin\\Dropbox\\WES CNV\\calling results 3 methods\\LDseq\\Choose reference\\summary")
write.csv(prop.table,file="prop.table.csv")

