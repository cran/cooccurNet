
#' @importFrom   foreach %do%
cooccur.dataprepreprocess.preprocess <- function(data=list(), threshold=0.9, parallel=FALSE, cpus=NA, memory=NA, debug=FALSE ){
  if(!requireNamespace("foreach", quietly = TRUE)){
    stop("Package 'foreach' is required.")
  }

  #check param 1
  nrow = nrow(data$matrix)
  ncol = ncol(data$matrix)


  if(nrow<=1 | ncol <=1){
    stop("The matrix should have more than one row and one column.")
  }


  #initialize frequency matrix of the original
  data$freqMatrix = matrix(nrow=nrow,ncol=ncol)
  storage.mode(data$freqMatrix) <- "integer"

  seqlevel = c()
  fable <- c()


  #sfInit(parallel=TRUE, cpus=cooccur.detectCores())
  #print(sprintf('%s cpus to be used', sfCpus()))

  #step1. calculate frequency
  t = Sys.time()
  #data$freqMatrix<-sfApply(data$matrix, 2, function(x) {x= as.character(x);table(x)[x]})
  #data$freqMatrix<-apply(data$matrix, 2, function(x) {x= as.character(x);table(x)[x]})
  #for(i in 1:ncol){
  i = 1
  foreach::foreach(i =  1:ncol) %do% {
    #print(i)
    x = as.character(data$matrix[,i])
    data$freqMatrix[,i] = table(x)[x]
    #rm(x)
  }
  cooccur.printTimeCost('getFrequence data time cost',t,debug)

  #step2. filter matrix columns by parameter "threshold"
  t = Sys.time()
  fable = c()
  storage.mode(fable) <- "integer"
  #fable<-sfApply(data$freqMatrix, 2, function(x){ifelse(round(max(x)/nrow,5)>=threshold, 1,0)})
  #fable<-apply(data$freqMatrix, 2, function(x){ifelse(round(max(x)/nrow,5)>=threshold, 1,0)})
  #for(i in 1:ncol){
  #print(threshold)
  i = 1
  foreach::foreach(i =  1:ncol) %do% {
    #print(i)
    x =  data$freqMatrix[,i]
    #print(round(max(x)/nrow,5))
    #print(threshold)
    if(round(max(x)/nrow,5)>=threshold){
      fable <- c(fable,i)
      #print(fable)
    }
  }
  #print(length(fable)/ncol(data$matrix))
  #print(length(fable))

  #2016-08-30 begin
  #if(is.na(memory)==TRUE){
    #print(length(fable)/ncol(data$matrix))
    #print(ncol(data$matrix))
    #if(ncol(data$matrix) <= 100){
      #memory="memory"
    #}else if(length(fable)/ncol(data$matrix)<=0.3){
      #memory="sparse"
    #}else{
      #memory="memory"
    #}
  #}
  #2016-08-30 end







  if(length(fable)>0){
    data$matrix <- data$matrix[,-fable]
    data$freqMatrix  <- data$freqMatrix[,-fable]
    data$original <- data$original[,-fable]
  }
  #sfStop()


  inmemoryflag = inmemoryflag(nrow(data$matrix), ncol(data$matrix))

  if(is.na(memory)==TRUE){
    if(inmemoryflag==TRUE){
      memory="memory"
    }else{
      memory="sparse"
    }
  }

  msg = paste('using memory,', memory, sep="")
  cooccur.printTimeCost(msg, t, debug)

  cooccur.printTimeCost('compareWiththreshold data time cost',t, debug)
  rm(fable)



  #step3. create big memory matrix of column-column relation
  t = Sys.time()
  filename = data$dt_idxtable_filename

  #2016-09-19 begin
  #data$dt_idxtable <- cooccur.geneColumnCombination.bigmemory(ncol(data$matrix),filename,memory)	#recalculate ncol(matrix)

  data$dt_idxtable <- cooccur.geneColumnCombination.normal(ncol(data$matrix),filename,memory)
  #2016-09-19 end
   #data$dt_idxtable <- cooccur.geneColumnCombination.dt(ncol(data$matrix))
  cooccur.printTimeCost('create ColumnCombination table dt_idxtable time cost',t,debug)

  #gc memory
  gc()

  #memory="memory"

  #step4. calculation bigramfrequency
  t = Sys.time()
  if(memory=="memory"){
    cat("calculating bigram frequency ......")
    message("")
    data$bigramFreqList <- cooccur.dataprepreprocess.bigramfrequence(data)
    #print(data$bigramFreqList[,])
    data$memory <- "memory"
  }else if(memory=="sparse"){
    data$memory <- "sparse"
  }



  cooccur.printTimeCost('create bi-grams frequence table time cost',t,debug)
  #print(paste('create bi-grams frequence table time cost:', round((Sys.time()-t),2),' secs')	)


  #return object
  return(data)


}


#' @importFrom   foreach %do%
#' @importFrom utils memory.size setTxtProgressBar txtProgressBar write.table
cooccur.dataprepreprocess.bigramfrequence <- function(sequences,colsperIter=20){

  #print(dim(sequences$matrix))
  #print(dim(sequences$original))


  if(!requireNamespace("bigmemory", quietly = TRUE)){
    stop("Package 'bigmemory' is required.")
  }
  if(!requireNamespace("foreach", quietly = TRUE)){
    stop("Package 'foreach' is required.")
  }
  #t1 = Sys.time()
  len = length(sequences$dt_idxtable[,2])
  nrow = nrow(sequences$matrix)

  bigramfreqx = bigmemory::big.matrix(nrow=len, ncol=nrow, init=NA, dimnames=list(NULL,1:nrow))
  #storage.mode(bigramfreqx) <- "integer"
  #cooccur.printTimeCost('cooccur.dataprepreprocess.bigramfrequence',t1,debug)

  #print(len)

  #2016-09-05
  colsperIter = 50
  #if(len<2001){
    #colsperIter = 200
  #}else if(len<100001){
    #colsperIter = 10000
  #}else if(len<1000001){
    #colsperIter = 20
  #}
  Iter = ceiling(len / colsperIter)
  start <- 1
  end <- 1
  idx <- 1
  pb <- txtProgressBar(style = 3)
  if(Iter == 1){
    progress = seq(1, 1)
  }else{
    progress = seq(0, 1, 1/(Iter-1))
  }

  #for(k in 1:Iter){
  k = 1
  foreach::foreach(k =  1:Iter) %do% {
    bigramfreq <- c()
    end <- k * colsperIter
    if(end > len) end = len

    i <- sequences$dt_idxtable[,2][start:end]
    j <- sequences$dt_idxtable[,3][start:end]



    #t = Sys.time()
    #print(dim(sequences$matrix[,i]))
    xx = matrix(nrow=nrow, ncol=length(i))
    storage.mode(xx) <- "integer"



    #2016-09-04
    #xx <- paste(sequences$original[,i], sequences$original[,j],sep = "")
    #xx =  big.matrix(nrow=nrow, ncol=length(i), init=NA, dimnames=list(NULL,1:length(i)))
    xx =   (round( sequences$matrix[,i] / sequences$matrix[,j],5) + sequences$matrix[,i]) * 100000
	  #xx =    sequences$matrix[,i] * sequences$matrix[,j] + sequences$matrix[,i]
    xx = apply(xx,2, function(x){sequences$constantList$biseqlevel[match(x, sequences$constantList$biseqidlevel)]})
    #sequences$constantList$biseqlevel[match(aa, sequences$constantList$biseqidlevel)]

    #cooccur.printTimeCost('(matrix(bigramSeqs,nrow=nrow))',t)

    #t = Sys.time()
    #aa = t(apply(xx, 2, function(x) {x= as.character(x);table(x)[x]}))
    aa = t(apply(xx, 2, function(x) {table(x)[x]}))
    #cooccur.printTimeCost('lapply(seq_len(ncol(xx)), function(x){ f[[x]][xx[,x]]})',t)

    bigramfreqx[start:end, ] <- aa

    start <- end + 1
    end <- k * colsperIter

    setTxtProgressBar(pb, progress[k])
  }
  close(pb)
  return(bigramfreqx)
}

###
