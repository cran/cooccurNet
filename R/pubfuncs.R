#' @importFrom parallel detectCores
cooccur.detectCores <- function(){
  #require(parallel)
	if(!requireNamespace("parallel", quietly = TRUE)){
		stop("Package 'parallel' is required.")
	}
	cpus <- parallel::detectCores()
	#sysname = as.character(Sys.info()["sysname"])
	#cpus = 1
	#if(sysname=="Windows"){
		#cpus <- detectCores(logical=F)
	#}else if(sysname=="Linux"){
		#cpus <- detectCores()
	#}
	#print(cpus)
	return(cpus)
}


cooccur.removeidxtablefile <- function(filename){
	#print(filename)
	#binfile = paste(filename, ".bin", sep="")
	descfile = paste(filename, ".desc", sep="")

	#if(file.exists(binfile)){
		#file.remove(binfile)
	#}
	if(file.exists(descfile)){
		#file.remove(descfile)
	}
}



cooccur.geneColumnCombination.bigmemory <- function(ncol,filename,memory){
	#require(bigmemory)
	if(!requireNamespace("bigmemory", quietly = TRUE)){
		stop("Package 'bigmemory' is required.")
	}
  gc()

  #2016-09-04
  if(memory=="sparsex"){
  #if(memory=="sparse"){
  	#print(filename)
  	binfile = paste(filename, ".bin", sep="")
  	descfile = paste(filename, ".desc", sep="")

  	if(file.exists(binfile)){
  		file.remove(binfile)
  	}
  	if(file.exists(descfile)){
  		file.remove(descfile)
  	}

  	nrow = ncol * (ncol - 1) / 2
  	dt_idxtable <- bigmemory::filebacked.big.matrix(nrow=nrow, ncol=3, type="integer", init=NA, dimnames=list(NULL,1:3), backingfile=binfile, descriptorfile=descfile)

  	dt_idxtable[,1] = seq(1,nrow)
  	dt_idxtable[,2] = as.vector(unlist(sapply(seq_len(ncol-1),function(x) a <- rep(x, times=(ncol-x)))))
  	dt_idxtable[,3] = as.vector(unlist(sapply(2:(ncol),function(x) a <- x : ncol )))
  	return(dt_idxtable)

  #2016-09-04
  #}else if(memory=="memory"){
  }else{
    nrow = ncol * (ncol - 1) / 2
    dt_idxtable <- bigmemory::big.matrix(nrow=nrow, ncol=3, type="integer", init=NA, dimnames=list(NULL,1:3))

    dt_idxtable[,1] = seq(1,nrow)
    dt_idxtable[,2] = as.vector(unlist(sapply(seq_len(ncol-1),function(x) a <- rep(x, times=(ncol-x)))))
    dt_idxtable[,3] = as.vector(unlist(sapply(2:(ncol),function(x) a <- x : ncol )))
    return(dt_idxtable)
  }
	#print(storage.mode(dt_idxtable))

}

cooccur.geneColumnCombination.bigmemory.old <- function(ncol,filename,memory){
  #require(bigmemory)
  if(!requireNamespace("bigmemory", quietly = TRUE)){
    stop("Package 'bigmemory' is required.")
  }
  gc()
  
  if(memory=="sparse"){
    #print(filename)
    binfile = paste(filename, ".bin", sep="")
    descfile = paste(filename, ".desc", sep="")
    
    if(file.exists(binfile)){
      file.remove(binfile)
    }
    if(file.exists(descfile)){
      file.remove(descfile)
    }
    
    nrow = ncol * (ncol - 1) / 2
    dt_idxtable <- bigmemory::filebacked.big.matrix(nrow=nrow, ncol=3, type="integer", init=NA, dimnames=list(NULL,1:3), backingfile=binfile, descriptorfile=descfile)
    
    dt_idxtable[,1] = seq(1,nrow)
    dt_idxtable[,2] = as.vector(unlist(sapply(seq_len(ncol-1),function(x) a <- rep(x, times=(ncol-x)))))
    dt_idxtable[,3] = as.vector(unlist(sapply(2:(ncol),function(x) a <- x : ncol )))
    return(dt_idxtable)
    
  }else if(memory=="memory"){
    nrow = ncol * (ncol - 1) / 2
    dt_idxtable <- bigmemory::big.matrix(nrow=nrow, ncol=3, type="integer", init=NA, dimnames=list(NULL,1:3))
    
    dt_idxtable[,1] = seq(1,nrow)
    dt_idxtable[,2] = as.vector(unlist(sapply(seq_len(ncol-1),function(x) a <- rep(x, times=(ncol-x)))))
    dt_idxtable[,3] = as.vector(unlist(sapply(2:(ncol),function(x) a <- x : ncol )))
    return(dt_idxtable)
  }
  #print(storage.mode(dt_idxtable))
  
}


cooccur.printTimeCost <- function(msg, t1, debug=FALSE){
  #print(t)
	if(debug==TRUE){
		message("")
		msg = paste("    ",msg,sep="")
		message(paste(msg, ":", round(difftime(Sys.time(),t1, units="secs"),2),' secs;memory applied:',memory.size(T), "Mb",sep=""))
	}

}

cooccur.writelist <- function(list,file){
	#write.list(list,filename = file)
	for(i in 1:length(list)){
		write.table(paste(unlist(list[i]),collapse =  "	"), file=file, append = TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
	}
}

cooccur.writeMatrix <- function(m,file){
	#write.list(list,filename = file)
	for(i in 1:nrow(m)){
		write.table(paste(m[i,],collapse =  "	"), file=file, append = TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
	}
}

cooccur.writetable <- function(df,file){
	write.table(df, file=file, append = TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
}

#' @importFrom   foreach %do%
cooccur.constant <- function(){
	if(!requireNamespace("foreach", quietly = TRUE)){
		stop("Package 'foreach' is required.")
	}
  seqlevel = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","-","0","1","2","3","4","5","6","7","8","9")
  seqidlevel = c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,73, 79, 83, 89, 97, 101, 103, 107, 109, 113,127,131,137,139,149,151,157,163,173,179,181,191,193  )

	biseqlevel = c()
	biseqidlevel = c()
	#for(i in 1:length(seqlevel)){
	i = 1
	foreach::foreach(i =  1:length(seqlevel)) %do% {
		#for(j in 1:length(seqlevel)){
		j = 1
		foreach::foreach(j =  1:length(seqlevel)) %do% {
			a = seqlevel[i]
			b = seqlevel[j]
			c = paste(a,b,sep="")
			a1 = seqidlevel[i]
			b1 = seqidlevel[j]
			c1 = (round( a1 / b1, 5) + a1 ) * 100000
			#c1 =  a1 * b1  + a1
			biseqlevel = c(biseqlevel, c)
			biseqidlevel = c(biseqidlevel, c1)
		}
	}
	return(list(biseqlevel=biseqlevel, biseqidlevel=biseqidlevel))
}


cooccur.conservativeFilter <- function(conservativeFilter,dataType){
  if(is.na(conservativeFilter)){
    if(dataType=="DNA"){
      conservativeFilter = 0.95
    }else if(dataType=="protein"){
      conservativeFilter = 0.95
    }else if(dataType=="SNP"){
      conservativeFilter = 0.95
    }else{
      conservativeFilter = 0.95
    }
  }else if(is.null(conservativeFilter)){
    if(dataType=="DNA"){
      conservativeFilter = 0.95
    }else if(dataType=="protein"){
      conservativeFilter = 0.95
    }else if(dataType=="SNP"){
      conservativeFilter = 0.95
    }else{
      conservativeFilter = 0.95
    }    
  }else{
    x = as.numeric(conservativeFilter)
    if(x<0 | x>1){
      stop("the param 'conservativeFilter' is invalid which should in the range of [0,1]")
    }
  }  
  #print(sprintf("conservativeFilter:%s",conservativeFilter))
  return(conservativeFilter)
  
}


cooccur.cooccurFilter <- function(cooccurFilter,dataType){
  #print(is.null(cooccurFilter))
  if(is.null(cooccurFilter)){
    if(dataType=="DNA"){
      cooccurFilter = 1
    }else if(dataType=="protein"){
      cooccurFilter = 0.9
    }else if(dataType=="SNP"){
      cooccurFilter = 1
    }else{
      cooccurFilter = 1
    }
  }else if(is.na(cooccurFilter)){
    if(dataType=="DNA"){
      cooccurFilter = 1
    }else if(dataType=="protein"){
      cooccurFilter = 0.9
    }else if(dataType=="SNP"){
      cooccurFilter = 1
    }else{
      cooccurFilter = 1
    }
  }else{
    x = as.numeric(cooccurFilter)
    if(x<0 | x>1){
      stop("the param 'cooccurFilter' is invalid which should in the range of [0,1]")
    }
  }  
  #print(sprintf("cooccurFilter:%s",cooccurFilter))
  return(cooccurFilter)
  
}

#' @importFrom utils memory.limit 
inmemoryflag <- function(nrow, ncol){
  alpha = 0.5
  memsize = memory.limit()
  pred = 12 * nrow * ncol * (ncol - 1) / 1024 / 1024
  #print(pred)
  #print(memsize *alpha)
  if(pred >= memsize *alpha){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

