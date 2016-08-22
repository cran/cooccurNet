
#' @importFrom   foreach %do%
cooccur.networkpvalue.calculateNetWorkPvalue <- function(cooccurrence, sequences, ptimes = 100,alpha=0.9,parallel=FALSE,cpus=cooccur.detectCores(),filterfile=NA,rawfile=NA,modulefile=NA, propertyfile=NA, cooccurfile=NA){
	if(!requireNamespace("Matrix", quietly = TRUE)){
		stop("Package 'Matrix' is required.")
	}
	if(!requireNamespace("foreach", quietly = TRUE)){
		stop("Package 'foreach' is required.")
	}
	pvalue = cooccurrence[,3]
	#print(pvalue)
	message("")

	steps = 0
	#for(i in 1:ptimes){
	i = 1
	foreach::foreach(i =  1:ptimes) %do% {
		steps = steps + 1
		cat(sprintf("    calculating shuffled network (%s) ......", steps))
		message("")
		sequences$matrix <- cooccur.networkpvalue.shuffleMatrix(sequences)
		if(sequences$memory == "memory"){
			sequences$bigramFreqList <- cooccur.dataprepreprocess.bigramfrequence(sequences,colsperIter=20)
		}
		cooccurr <- cooccur.networkpvalue.networkcooccurrence(sequences,alpha,parallel,cpus)
		#print(dim(cooccurr))
		#cooccurrence = cbind(cooccurrence,cooccurr)
		cooccurrence = Matrix::cBind(cooccurrence,cooccurr)
	}
	#print(cooccurrence)
	#print(cooccurrence[,3:ncol(cooccurrence)])
	cooccurr = t(apply(cooccurrence[,3:ncol(cooccurrence)],1,sort,decreasing=T))
	pvalues = sapply(seq_len(length(pvalue)),cooccur.networkpvalue.pvalueCompare,pvalue, cooccurr )
	#print(pvalues)
	#print(unlist(pvalues))
	return(data.frame(i=cooccurrence[,1],j=cooccurrence[,2],cooccur=cooccurrence[,3],pvalue=pvalues))
	#print(cooccurr)

}




cooccur.networkpvalue.networkcooccurrence <- function(sequences,alpha=0.9,parallel=FALSE,cpus=cooccur.detectCores()){
	cooccurrence = c()
	df_cooccurrence = cooccur.gennetework.cooccurnetworks(sequences,steps=NA)
	colnames(df_cooccurrence) =seq(1:ncol(df_cooccurrence))
	rownames(df_cooccurrence) = seq(1:nrow(sequences$matrix))
	#cooccurrence <- lapply(seq_len(ncol(df_cooccurrence)), cooccur.networkpvalue.getcooccurrence, sequences, df_cooccurrence)
	#cooccurrence = do.call("rbind",cooccurrence)

	cooccurrence <- cooccur.gennetework.outputNetWorkcooccurrence(sequences, df_cooccurrence, shuffeld=TRUE)

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


cooccur.networkpvalue.shuffleMatrix <- function(sequences){
	nrow = nrow(sequences$matrix)
	smatrix = apply(sequences$matrix,2, function(x){
		nrow = length(x)
		x = x[sample(seq_len(nrow),nrow,replace=FALSE)]
	})
	#print(smatrix)
	return(smatrix)
}

cooccur.networkpvalue.pvalueCompare <- function(x, pvalue, cooccurr){
	  a = pvalue[x]
	  b = cooccurr[x,]
	  pvalue = round((min(which(a==b))-1)/(ncol(cooccurr)-1),3)
	  pvalue = ifelse(pvalue==0,"<0.01",pvalue)
	  return(pvalue)
}


