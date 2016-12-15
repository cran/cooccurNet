
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


