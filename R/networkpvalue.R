
#' @importFrom   foreach %do%
cooccur.networkpvalue.calculateNetWorkPvalue <- function(df_cooccurrence, cooccurrence=NA, sequences, pvaluefile=NA, networkpfile=NA,ptimes = 100,alpha=0.9, parallel=FALSE, debug=FALSE){
#cooccurrence, sequences, ptimes, alpha, debug
  pvalues = list(siteco=NULL, networks=NULL)
  if(!requireNamespace("Matrix", quietly = TRUE)){
		stop("Package 'Matrix' is required.")
	}
	if(!requireNamespace("foreach", quietly = TRUE)){
		stop("Package 'foreach' is required.")
	}
  #if(is.na(df_cooccurrence) || is.null(df_cooccurrence)){
    #stop("'df_cooccurrence' cannot be NA.")
  #}
  #######################################################


  networkdf_cooccurrence = Matrix::Matrix(nrow=nrow(df_cooccurrence),ncol=ncol(df_cooccurrence), sparse=TRUE, data=0)
  #print(object.size(networkdf_cooccurrence))
  #print(dim(networkdf_cooccurrence))

  binomflag = FALSE
  if(nrow(df_cooccurrence)<=100){
    #binomflag = FALSE
  }else{
    #binomflag = TRUE
  }

  if(binomflag==TRUE && !is.na(pvaluefile)){
    pvalues$siteco = (data.frame(i=cooccurrence[,1],j=cooccurrence[,2],cooccur=cooccurrence[,3],pvalue=coocur.gennetwork.binom.test(df_cooccurrence)))
    if(is.na(networkpfile)){
      return(pvalues)
    }
  }
  #
	#pvalue = cooccurrence[,3]
	#print(pvalue)
	message("")
	#print(str(sequences))

	steps = 0
	for(i in 1:ptimes){
	#i = 1
	#foreach::foreach(i =  1:ptimes) %do% {
		steps = steps + 1
		cat(sprintf("    calculating shuffled network (%s) ......", steps))
		message("")
		t = as.character(Sys.time())
		sequences$matrix = cooccur.networkpvalue.shuffleMatrix(sequences,debug)

		cooccur.printTimeCost("cooccur.networkpvalue.shuffleMatrix ",t,debug)

		if(sequences$memory == "memory"){
			sequences$bigramFreqList = cooccur.dataprepreprocess.bigramfrequence(sequences, colsperIter=20, debug)
		}

		t = as.character(Sys.time())


		#2016-11-23 begin
		#cooccurr = cooccur.networkpvalue.networkcooccurrence(sequences,alpha, parallel, debug)


		networkcooccurr = c()
		#sequences, alpha=0.9, steps, debug=FALSE
		#2016-09-26 begin
		df_cooccurrence = cooccur.gennetework.cooccurnetworks(sequences, alpha, steps=NA, parallel, debug)
		#df_cooccurrence = cooccur.networkpvalue.cooccurnetworks(sequences, alpha, steps=NA, debug)
		#2016-09-26 end

    if(!is.na(networkpfile) ){
      networkdf_cooccurrence = networkdf_cooccurrence + df_cooccurrence
    }

		if(!is.na(pvaluefile) && binomflag==FALSE){
		  colnames(df_cooccurrence) =seq(1:ncol(df_cooccurrence))
		  rownames(df_cooccurrence) = seq(1:nrow(sequences$matrix))
		  #cooccurrence <- lapply(seq_len(ncol(df_cooccurrence)), cooccur.networkpvalue.getcooccurrence, sequences, df_cooccurrence)
		  #cooccurrence = do.call("rbind",cooccurrence)
		  t = Sys.time()
		  networkcooccurr = cooccur.gennetework.outputNetWorkcooccurrence(sequences, df_cooccurrence, shuffeld=TRUE, parallel = parallel, debug=debug)
		  #print(cooccurrence)
		  cooccur.printTimeCost('cooccur.gennetework.outputNetWorkcooccurrence cooccur time cost',t,debug)

		  #2016-11-23 end

		  cooccur.printTimeCost("cooccur.networkpvalue.networkcooccurrence ",t,debug)

		  sequences$bigramFreqList = NULL
		  #print(dim(cooccurr))
		  #cooccurrence = cbind(cooccurrence,cooccurr)
		  #print(cooccurrence)
		  cooccurrence = Matrix::cBind(cooccurrence,networkcooccurr[,1])
		  rm(networkcooccurr)
		}
    rm(df_cooccurrence)

		gc()
	}
	##print(cooccurrence)
	##print(cooccurrence[,3:ncol(cooccurrence)])
	if(!is.na(pvaluefile) && binomflag==FALSE){
	  cooccurr = t(apply(cooccurrence[,3:ncol(cooccurrence)],1,sort,decreasing=T))
	  pvalues$siteco = (data.frame(i=cooccurrence[,1],j=cooccurrence[,2],cooccur=cooccurrence[,3],pvalue = sapply(seq_len(length(cooccurrence[,3])),cooccur.networkpvalue.pvalueCompare, cooccurrence[,3], cooccurr )))

	}
	if(!is.na(networkpfile) ){
	  #print("---------------------------------------------------------")
	  #print(df_cooccurrence)
	  #print("--------------------------------------------------------")
	  pvalues$networks =  (networkdf_cooccurrence / ptimes)
	  #print(pvalues$networks)
	}

	rm(networkdf_cooccurrence)
	rm(cooccurrence)
	#rm(df_cooccurrence)
	rm(sequences)
	gc()
	return(pvalues)
	###print(pvalue)
	##pvalues = sapply(seq_len(length(pvalue)),cooccur.networkpvalue.pvalueCompare,pvalue, cooccurr )
	##pvalues = sapply(seq_len(length(pvalue)),cooccur.networkpvalue.pvalueCompare,pvalue, cooccurr )
	###print(pvalues)
	###print(unlist(pvalues))
	##return(data.frame(i=cooccurrence[,1],j=cooccurrence[,2],cooccur=cooccurrence[,3],pvalue=pvalues))
	#return(data.frame(i=cooccurrence[,1],j=cooccurrence[,2],cooccur=cooccurrence[,3],pvalue = sapply(seq_len(length(pvalue)),cooccur.networkpvalue.pvalueCompare,pvalue, cooccurr )))
  ###print(cooccurr)

}

#' @importFrom   foreach %do%
#' @importFrom   foreach %dopar%
cooccur.networkpvalue.calculateNetWorkPvalue.new <- function(df_cooccurrence=NA, cooccurValues=NA, sequences, pvaluefile=NA, networkpfile=NA,ptimes = 100,alpha=0.9, parallel=FALSE, debug=FALSE){
  #cooccurrence, sequences, ptimes, alpha, debug
  pvalues = list(siteco=NULL, networks=NULL)
  if(!requireNamespace("Matrix", quietly = TRUE)){
    stop("Package 'Matrix' is required.")
  }
  if(!requireNamespace("foreach", quietly = TRUE)){
    stop("Package 'foreach' is required.")
  }

  #######################################################
  #require(doParallel)

  #networkdf_cooccurrence = Matrix::Matrix(nrow=nrow(df_cooccurrence),ncol=ncol(df_cooccurrence), sparse=TRUE, data=0)


  cores <- cooccur.detectCores()
  ppmatrix = splitArrayByCores(ptimes, cores$cpus)
  if(nrow(ppmatrix)<cores$cpus){
    cores$cpus = nrow(ppmatrix)
  }

  if(parallel==TRUE){
    #cores = 3
    #cl <- makeCluster(cores)
    #registerDoParallel(cl)
    #if(cores$win==TRUE){
      #print("cl = parallel::makeCluster(cores$cpus)")
      #cl = parallel::makeCluster(cores$cpus)
    #}else{
      #print(paste("cores in use:",cores$cpus,sep=""))
      #print("cl = parallel::makeCluster(cores$cpus, type = 'FORK')")
      #cl = parallel::makeCluster(cores$cpus, type = "FORK")
    #}

    print(paste("cores in use:",cores$cpus,sep=""))
    print("cl = parallel::makeCluster(cores$cpus)")
    cl = parallel::makeCluster(cores$cpus)
    doParallel::registerDoParallel(cl)
  }

  #2017-01-11 begin
  #df_cooccurrence_vector = cooccur.networkpvalue.dfcooccurrencevector(sequences$matrix,sequences,alpha)
  #len = nrow(sequences$dt_idxtable)
  #cooccurValues = matrix(nrow=len,ncol=0)
  #starts = colnames(sequences$matrix)[sequences$dt_idxtable[,2]]
  #ends = colnames(sequences$matrix)[sequences$dt_idxtable[,3]]
  #cooccurValues = cbind(cooccurValues, starts)
  #cooccurValues = cbind(cooccurValues, ends)
  #cooccurValues = cbind(cooccurValues, df_cooccurrence_vector)
  #2017-01-11 end

  #tryCatch({
    i = 0
    if(parallel==TRUE){
      #library(doParallel)
      ptimescooccurValues = foreach::foreach(i =  1:cores$cpus, .combine=cbind) %dopar% {
        x = ppmatrix[i,]
        return(doParallelPtimes(x,sequences,debug, alpha))
      }
      parallel::stopCluster(cl)
    }else{
      ptimescooccurValues = foreach::foreach(i =  1:cores$cpus, .combine=cbind) %do% {
        x = ppmatrix[i,]
        return(doParallelPtimes(x,sequences,debug, alpha))
      }
    }


    cooccurValues = cbind(cooccurValues, ptimescooccurValues)


    #print(cooccurValues)
    #print(dim(cooccurValues))

    cooccurr = t(apply(cooccurValues[,3:ncol(cooccurValues)],1,sort,decreasing=T))
    pvalues$siteco = (data.frame(i=cooccurValues[,1],j=cooccurValues[,2],cooccur=cooccurValues[,3],pvalue = sapply(seq_len(length(cooccurValues[,3])),cooccur.networkpvalue.pvalueCompare, cooccurValues[,3], cooccurr )))

    rm(cooccurValues)
    rm(ptimescooccurValues)
    rm(sequences)
    gc()
    return(pvalues)

  #}, error=function(e){
      #parallel::stopCluster(cl)
      #print(e)
    #}
  #)

  gc()

}


doParallelPtimes <- function(x, sequences,debug, alpha){
  message("")
  print(x)
  len = nrow(sequences$dt_idxtable)
  cooccurValues = matrix(nrow=len,ncol=0)
  for(idx in x[1]:x[2]){
    shuffleMatrix = cooccur.networkpvalue.shuffleMatrix(sequences,debug)


    df_cooccurrence_vector = cooccur.networkpvalue.dfcooccurrencevector(shuffleMatrix,sequences,alpha,debug)

    cooccurValues = cbind(cooccurValues,df_cooccurrence_vector)
    #comparedfile = "D:\\foreach_test\\"
    #xxxfile = paste(comparedfile,idx,".txt",sep="")
    #if(file.exists(xxxfile)){
      #file.remove(xxxfile)
    #}
    #write.table(df_cooccurrence_vector, file=xxxfile, append = TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
  return(cooccurValues)
}

cooccur.networkpvalue.dfcooccurrencevector <- function(shuffleMatrix,sequences,alpha,debug){
  t = as.character(Sys.time())

  len = nrow(sequences$dt_idxtable)
  df_cooccurrence_vector = c()
  length(df_cooccurrence_vector) <- len

  nrow = nrow(shuffleMatrix)

  colsperIter = 5
  Iter = ceiling(len / colsperIter)
  start <- 1
  end <- 1

  k = 1
  for(k in 1:Iter){
    end <- k * colsperIter
    if(end > len) end = len
    i = sequences$dt_idxtable[,2][start:end]
    j = sequences$dt_idxtable[,3][start:end]
    x = round(sequences$freqMatrix[,i]/nrow,5)
    y = round(sequences$freqMatrix[,j]/nrow,5)
    xMy = x * y
    xy = sqrt(alpha * xMy)

    a = shuffleMatrix[,i]
    b = shuffleMatrix[,j]
    xx =   (round( a / b, 5) + a ) * 100000

    if(is.vector(xx)==FALSE){
      xx = apply(xx,2, function(x){sequences$constantList$biseqlevel[match(x, sequences$constantList$biseqidlevel)]})
      xx = (apply(xx, 2, function(x) {table(x)[x]}))
    }else{
      xx = sequences$constantList$biseqlevel[match(xx, sequences$constantList$biseqidlevel)]
      xx = table(xx)[xx]
      xx = as.matrix(xx)
    }

    xx = round(xx/nrow,5)
    xx = ifelse(xx>=xy,1,0)

    #df_cooccurrence =  Matrix::cBind(df_cooccurrence, xx)
    df_cooccurrence_vector = c(df_cooccurrence_vector, apply(xx,2, function(x){return(round(sum(x)/length(x),3))}))

    #print(df_cooccurrence_vector)
    rm(xx)
    rm(xy)

    start <- end + 1
    end <- k * colsperIter
  }

  cooccur.printTimeCost('cooccur.networkpvalue.dfcooccurrencevector time cost',t,debug)
  return(df_cooccurrence_vector)
}


#sequences,alpha,debug
cooccur.networkpvalue.networkcooccurrence <- function(sequences,alpha=0.9, parallel=FALSE, debug=FALSE){
	cooccurrence = c()
	#sequences, alpha=0.9, steps, debug=FALSE
	#2016-09-26 begin
	df_cooccurrence = cooccur.gennetework.cooccurnetworks(sequences, alpha, steps=NA, parallel, debug)
	#df_cooccurrence = cooccur.networkpvalue.cooccurnetworks(sequences, alpha, steps=NA, debug)
	#2016-09-26 end
	colnames(df_cooccurrence) =seq(1:ncol(df_cooccurrence))
	rownames(df_cooccurrence) = seq(1:nrow(sequences$matrix))
	#cooccurrence <- lapply(seq_len(ncol(df_cooccurrence)), cooccur.networkpvalue.getcooccurrence, sequences, df_cooccurrence)
	#cooccurrence = do.call("rbind",cooccurrence)
	t = Sys.time()
	cooccurrence = cooccur.gennetework.outputNetWorkcooccurrence(sequences, df_cooccurrence, shuffeld=TRUE, parallel = parallel, debug=debug)
	#print(cooccurrence)
	cooccur.printTimeCost('cooccur.gennetework.outputNetWorkcooccurrence cooccur time cost',t,debug)

	return(cooccurrence)
}


cooccur.networkpvalue.getcooccurrence.old <- function(colid,sequences,df_cooccurrence){
	x = df_cooccurrence[,colid]
	cooccurrence = 0
	len = length(which(x==1))
	sub_dfcooc <- subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]==colid)
	start <- colnames(sequences$matrix)[sub_dfcooc[2]]
	end <- colnames(sequences$matrix)[sub_dfcooc[3]]

	if(len>0){
		cooccurrence = round(len/length(x),3)
	}else{
		cooccurrence = 0.000
	}
	return(cooccurrence)
}

#' @importFrom   foreach %do%
cooccur.networkpvalue.shuffleMatrix <- function(sequences,debug=FALSE){
  t = as.character(Sys.time())
  smatrix = matrix(nrow=nrow(sequences$matrix), ncol=ncol(sequences$matrix), 0)

  #foreach::foreach(colid =  1:ncol(sequences$matrix)) %do% {
  for(colid in 1:ncol(sequences$matrix)){
    nrow = length(sequences$matrix[,colid])
    smatrix[,colid] = sequences$matrix[,colid][sample(seq_len(nrow),nrow,replace=FALSE)]
  }

  storage.mode(smatrix) = "integer"
  cooccur.printTimeCost('cooccur.networkpvalue.shuffleMatrix time cost',t,debug)
  rm(sequences)
  return(smatrix)
}


cooccur.networkpvalue.shuffleMatrix.old <- function(sequences,debug=FALSE){
  t = as.character(Sys.time())
	nrow = nrow(sequences$matrix)
	smatrix = apply(sequences$matrix,2, function(x){
		nrow = length(x)
		x = x[sample(seq_len(nrow),nrow,replace=FALSE)]
	})
	#print(smatrix)
	storage.mode(smatrix) = "integer"
	cooccur.printTimeCost('cooccur.networkpvalue.shuffleMatrix time cost',t,debug)
	return(smatrix)
}

cooccur.networkpvalue.pvalueCompare <- function(x, pvalue, cooccurr){
	  a = pvalue[x]
	  b = cooccurr[x,]
	  #print(a)
	  #print(b)
	  #print((min(which(b==a))-1))
	  #2016-09-20 begin
	  #pvalue = round((min(which(b==a))-1)/(ncol(cooccurr)-1),3)
	  #pvalue = round((min(which(a==b))-1)/(ncol(cooccurr)-1),3)
	  pvalue = round((max(which(b>=a))-1)/(ncol(cooccurr)-1),3)
	  #2016-09-20 end

	  #2016-12-13 begin
	  #pvalue = ifelse(pvalue==0,"<0.01",pvalue)

	  pvalue = ifelse(pvalue==0, paste("<", round(1/(ncol(cooccurr)-1),3), sep="") ,pvalue)
	  #2016-12-13 end

	  return(pvalue)
}


