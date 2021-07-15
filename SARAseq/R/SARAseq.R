SARAseq <-function(lrr, map, h1 = 5, h2 = 10, h3 = 15, alpha = 0.01, L = 100, type=NULL, outname){
  cplist <- vector("list",dim(lrr)[2])
  for (i in 1:dim(lrr)[2]) {
    modSaRa   = modifiedSaRa(lrr[,i], alpha = alpha, L = L, h1 = h1, h2 = h2, h3 = h3, type=typ)
    cplist[[i]]    = modSaRa$newcp
    cnv.state = modSaRa$cnv.state
    cnv.start = modSaRa$cnv.start
    cnv.end   = modSaRa$cnv.end
    cnv.n  <- length(cnv.state)
    if (length(cnv.n)==0) break
    output <- matrix(NA, cnv.n, 9)
    for (j in 1 : cnv.n) {
      output[j,1] <- i  #Save individual id
      index.start <- cnv.start[j]
      index.end   <- cnv.end[j]
      output[j,2] <- map[index.start,2] #chromosome number
      output[j,3] <- as.character(map[index.start,1]) #cnv start marker name
      output[j,5] <- as.numeric(map[index.start,3]) #cnv ending marker name
      output[j,4] <- as.character(map[index.end,1]) # cnv start position
      output[j,6] <- as.numeric(map[index.end,3]) #cnv ending position
      output[j,7] <- as.numeric(output[j,6]) - as.numeric(output[j, 5]) #length of CNV
      output[j,8] <- index.end - index.start + 1 #NumSNP
      output[j,9] <- as.character(cnv.state[j]) #copy number state, duplication or deletion
    }
    write.table(output, paste(outname, ".csv", sep=""), append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

  }
  return (cp = cplist)
}
