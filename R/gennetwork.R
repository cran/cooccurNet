cooccur.gennetework.calculateNetWork <- function(sequences=list(),alpha=0.9,parallel=FALSE, filterfile=NA,rawfile=NA,modulefile=NA, propertyfile=NA, cooccurfile=NA, pvaluefile=NA, networkpfile=NA, ptimes=100,debug=FALSE){

	steps = 0

	#cooccurrenceList = c()
	#cooccurrenceEdges = c()
	#cat(sprintf("alpha = %s", alpha))


	if(!is.na(filterfile)){
		if(file.exists(filterfile)){
			file.remove(filterfile)
		}
		t = Sys.time()
		steps = steps+1
		cat(sprintf("%s. creating filter file (%s) ......", steps,filterfile))
		t = Sys.time()
		cooccur.writetable(sequences$original,filterfile)
		cooccur.printTimeCost('write networks filterfile time cost',t,debug)
		cat("completed")
	}

	if(!is.na(rawfile)){
	  sequences$networkFile = rawfile
		if(file.exists(rawfile)){
			file.remove(rawfile)
		}
	}

	if(!is.na(networkpfile)){
	  #sequences$networkpfile = networkpfile
	  if(file.exists(networkpfile)){
	    file.remove(networkpfile)
	  }
	}

	if(!is.na(modulefile)){
	  sequences$moduleFile = modulefile
		if(file.exists(modulefile)){
			file.remove(modulefile)
		}
	}
	if(!is.na(cooccurfile)){
	  #sequences$siteCoFile = cooccurfile
		if(file.exists(cooccurfile)){
			file.remove(cooccurfile)
		}
		header = list()
		header$h = "Site_i	Site_j	Cooccur"
		cooccur.writetable(header$h,cooccurfile)
	}

	if(!is.na(propertyfile)){
	  sequences$propertyFile = propertyfile
		if(file.exists(propertyfile)){
			file.remove(propertyfile)
		}
		header = list()
		header$h = "name  Connectivity  Diameter  Radius  ConnectionEffcient"
		cooccur.writetable(header$h,propertyfile)
	}


	if(!is.na(pvaluefile)){
		sequences$siteCoFile = pvaluefile
		if(file.exists(pvaluefile)){
			file.remove(pvaluefile)
		}
		header = list()
		header$h = "Site_i Site_j Co-occur p-value"
		cooccur.writetable(header$h,pvaluefile)
	}


	message("")
	t = Sys.time()

	#2016-11-16 begin
	#df_cooccurrence = cooccur.gennetework.cooccurnetworks(sequences,alpha,steps,debug)
	df_cooccurrence = cooccur.gennetework.cooccurnetworks(sequences, alpha, steps, parallel, debug)

	#2016-11-16 end

	cooccur.printTimeCost('cooccur.gennetework.cooccurnetworks df_cooccurrence time cost',t, debug)
	#2016-09-19
	#print((df_cooccurrence))

	#2016-08-31
	#if(sequences$memory=="sparse") {
		steps = steps + 1
	#}
	colnames(df_cooccurrence) =seq(1:ncol(df_cooccurrence))
	rownames(df_cooccurrence) = seq(1:nrow(sequences$matrix))

  #print(paste(rawfile,"-rawfile",sep="--"))
	if(!is.na(rawfile)){
		steps = steps+1
		cat(sprintf("%s. creating network file (%s) ......", steps,rawfile))
		t = Sys.time()
		#cooccurrenceEdgesBM <- data.frame()

		#2016-09-21 begin
		#write rawfile comments only for test reason,
		lapply(seq_len(nrow(df_cooccurrence)), cooccur.gennetework.outputNetWork, sequences, df_cooccurrence, rawfile)


		cooccur.printTimeCost('write networks networkFile time cost',t, debug)
		#rm(cooccurrenceEdgesBM)
		gc()
		cat("completed")
	}
	message("")
	if(!is.na(modulefile)){
		steps = steps+1
		cat(sprintf("%s. creating module file (%s) ......", steps,modulefile))
		t = Sys.time()
		cooccurrenceModule <- lapply(seq_len(nrow(df_cooccurrence)), cooccur.gennetework.outputNetWorkModule,sequences, df_cooccurrence)
		cooccur.printTimeCost('create cooccurrence Module time cost',t,debug)
		t = Sys.time()
		cooccur.writelist(cooccurrenceModule,modulefile)
		cooccur.printTimeCost('write networks modulefile time cost',t,debug)
		cat("completed")
	}
	message("")
	if(!is.na(propertyfile)){
		steps = steps+1
		cat(sprintf("%s. creating property file (%s) ......",steps, propertyfile))
		t = Sys.time()
		cooccurrenceProperties <- lapply(seq_len(nrow(df_cooccurrence)), cooccur.gennetework.outputNetWorkProperties, sequences, df_cooccurrence)
		cooccur.printTimeCost('create cooccurrence property time cost',t, debug)
		t = Sys.time()
		cooccur.writelist(cooccurrenceProperties,propertyfile)
		cooccur.printTimeCost('write networks propertyfile time cost',t,debug)
		cat("completed")
	}
	message("")

	if(!is.na(cooccurfile)){
		cooccurrence = list()
		steps = steps+1
		cat(sprintf("%s. creating siteCoFile file (%s) ......",steps, cooccurfile))
		t = Sys.time()
		#cooccurrence <- lapply(seq_len(ncol(df_cooccurrence)), cooccur.gennetework.outputNetWorkcooccurrence, sequences, df_cooccurrence)
		cooccurrence <- cooccur.gennetework.outputNetWorkcooccurrence(sequences, df_cooccurrence, shuffeld=FALSE, parallel = parallel, debug=debug)
		#print(cooccurrence)
		cooccur.printTimeCost('cooccur.gennetework.outputNetWorkcooccurrence cooccurrence time cost',t,debug)
		t = Sys.time()
		#cooccur.writelist(cooccurrence,cooccurfile)
		cooccur.writeMatrix(cooccurrence,"cooccurfile")
		cooccur.printTimeCost('write networks siteCoFile time cost',t,debug)
		cat("completed")
	}
	message("")

	#2016-11-21 add networkpvalue
	#if pvaluefile or networkpvalue is not NA,

	if(!is.na(pvaluefile) || !is.na(networkpfile)){
		steps = steps+1



		#if(nrow(cooccurrence)==0){
		#if(length(cooccurrence)==0){
			#cooccurrence <- lapply(seq_len(ncol(df_cooccurrence)), cooccur.gennetework.outputNetWorkcooccurrence, sequences, df_cooccurrence)
			#print("-----------------------------------------------------------------------------------------")
			#sequences,df_cooccurrence,shuffeld=FALSE, debug
		if(!is.na(pvaluefile)){
		  cat(sprintf("calculating cooccurrence.........."))
		  message("")
		  t = Sys.time()
		  cooccurrence <- cooccur.gennetework.outputNetWorkcooccurrence(sequences, df_cooccurrence, shuffeld=FALSE, parallel = parallel, debug=debug)
		  colnames(cooccurrence) = c("Site_i","Site_j","Cooccur")
		  cooccur.printTimeCost('cooccur.gennetework.outputNetWorkcooccurrence cooccur time cost',t,debug)
		}else{
		  cooccurrence = NA
		}
		  #cooccurrence <- Matrix(cooccurrence)
		  #print(cooccurrence)


		  #rm(df_cooccurrence)
		  #gc()

		  #print(dim(cooccurrence))
		  #2016-09-21 begin
		  #print(object.size(cooccurrence)/1024/1024/1024)

		  #print(cooccurrence)
		#}#else{
			#cooccurrence = do.call("rbind", cooccurrence)
		#}


		t = Sys.time()
		#cooccurrence, sequences, ptimes = 100,alpha=0.9, debug=FALSE

		#2016-11-23 begin
		#binomflag = FALSE
		#if(nrow(df_cooccurrence)<=100){
		  ##binomflag = FALSE
		#}else{
		  ##binomflag = TRUE
		#}

		#if(binomflag==TRUE){
		  #pvalues = coocur.gennetwork.binom.test(df_cooccurrence)
		  ##print(length(pvalues))
		  ##print(length(cooccurrence[,1]))
		  ##print(length(cooccurrence[,2]))
		  ##print(length(cooccurrence[,3]))
		  #pvalues = (data.frame(i=cooccurrence[,1],j=cooccurrence[,2],cooccur=cooccurrence[,3],pvalue=pvalues))
		  ##pvalues = Matrix::cBind(cooccurrence, pvalues)

		#}else{
		if(!is.na(pvaluefile) && !is.na(networkpfile)){
		  cat(sprintf("%s. creating siteCoFile (%s) and networkEvaluate file (%s), sampleTimes: %s ......",steps, pvaluefile, networkpfile, ptimes))
		}else{
		  if(!is.na(pvaluefile)){
		    cat(sprintf("%s. creating siteCoFile (%s), sampleTimes: %s ......",steps, pvaluefile,ptimes))
		  }
		  if(!is.na(networkpfile)){
		    cat(sprintf("%s. creating networkEvaluate file (%s), sampleTimes: %s ......",steps, networkpfile,ptimes))
		  }
		}

		#cooccur.printTimeCost('before cooccur.networkpvalue.calculateNetWorkPvalue cooccur time cost',t,debug)
		pvalues = cooccur.networkpvalue.calculateNetWorkPvalue(df_cooccurrence, cooccurrence, sequences, pvaluefile=pvaluefile, networkpfile=networkpfile, ptimes=ptimes, alpha=alpha, parallel=parallel, debug=debug)
		#print(pvalues)
		#}

		#print(";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;")
		#print(df_cooccurrence)
		#print(";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;")

		#2016-10-18 begin
		#pvalues = cooccur.networkpvalue.calculateNetWorkPvalue(cooccurrence, sequences, ptimes, alpha, debug)
    #print(cooccurrence)
    #pvalues = coocur.gennetwork.binom.test(df_cooccurrence)
    #pvalues = (data.frame(i=cooccurrence[,1],j=cooccurrence[,2],cooccur=cooccurrence[,3],pvalue=pvalues))
    #coocur.gennetwork.binom.test(df_cooccurrence)
		#2016-10-18 end
		cooccur.printTimeCost('create cooccurrence p-value time cost',t,debug)
		#print(pvalues$siteco)
		if(!is.null(pvalues$siteco)){
		  t = Sys.time()
		  cooccur.writetable(pvalues$siteco,pvaluefile)
		  cooccur.printTimeCost('write networks siteCoFile time cost',t,debug)
		}

		#print(pvalues$networks)
		if(!is.null(pvalues$networks)){
		  t = Sys.time()
		  pvalues$networks = df_cooccurrence * pvalues$networks
		  lapply(seq_len(nrow(df_cooccurrence)), cooccur.gennetework.outputNetWorkpvalue, sequences, pvalues$network, networkpfile)
		  cooccur.printTimeCost('write networks networkpfile time cost',t,debug)

		  t = Sys.time()
		  newnetworkfile = paste(rawfile, "_new", sep="")
		  pvalues$networks = df_cooccurrence * pvalues$networks
		  lapply(seq_len(nrow(df_cooccurrence)), cooccur.gennetework.outputNewNetWorkFile, sequences, pvalues$network, newnetworkfile)
		  cooccur.printTimeCost('write networks newnetworkfile time cost',t,debug)

		}
		cat("completed")
		rm(cooccurrence)
		rm(df_cooccurrence)
		rm(pvalues)
		gc()
	}

	return(sequences)

}


#' @importFrom   foreach %do%
cooccur.gennetework.cooccurnetworks <- function(sequences, alpha=0.9, steps=0, parallel=FALSE, debug=FALSE){
  #print(alpha)

  #parallel = TRUE

	if(!requireNamespace("foreach", quietly = TRUE)){
		stop("Package 'foreach' is required.")
	}

  if(!requireNamespace("parallel", quietly = TRUE)){
    stop("Package 'parallel' is required.")
  }


  cores <- cooccur.detectCores()
  if(cores$win==TRUE){
    parallel = FALSE
  }

  df_cooccurrence = c()

	t = as.character(Sys.time())
	if(sequences$memory=="memory"){
	  #2016-08-30
	  #print(alpha)
	  if(!is.na(steps)){
	    steps = steps+1
	    cat(sprintf("%s. calculating networks ......", steps))
	    message("")
	  }

	  pb <- c()

	  if(debug){
	    pb <- txtProgressBar(style = 3)
	    progress = seq(0,1, 1/(nrow(sequences$dt_idxtable)-1))
	  }


		df_cooccurrence = data.frame()

		#2016-11-16 begin
		cooccurrenceList = c()
		#cooccurrenceList <- c(NA)
		length(cooccurrenceList) = nrow(sequences$dt_idxtable)
		#2016-11-16 end

		for(i in 1:nrow(sequences$dt_idxtable)){
			rowx = sequences$dt_idxtable[i,]

			cooccurrence = cooccur.gennetework.calucation(rowx,alpha,sequences)

			#print(cooccurrence)

			cooccurrenceList = c(cooccurrenceList,list(cooccurrence))


			if(debug){
  			#2016-08-30
  			setTxtProgressBar(pb, progress[i])
			}
		}
		if(debug){
  		#2016-08-30
  		close(pb)
		}
		#print(cooccurrenceList)
		df_cooccurrence = t(do.call(rbind, cooccurrenceList))
		rm(cooccurrenceList)
		gc()
		#cooccur.printTimeCost('cooccur.gennetework.cooccurnetworks time cost',t,debug)
		return(df_cooccurrence)
	#2016-09-19
	}else if(sequences$memory=="sparse" && parallel==TRUE){
		nrow = nrow(sequences$matrix)
		#df_cooccurrence = sparseMatrix(nrow, nrow(sequences$dt_idxtable), x=0)
		#bbbb = Matrix::Matrix(nrow=nrow,ncol=0,sparse=TRUE)
		#bbbb = list()
		if(!is.na(steps)){
			steps = steps+1
			cat(sprintf("%s. calculating networks ......", steps))
			message("")
		}

		cooccur.printTimeCost('begin parLapply ,,,,',t,debug)
		gc()
		M = 1:nrow(sequences$dt_idxtable)

		cores <- cooccur.detectCores()
		if(cores$win==TRUE){
		  print("cl = parallel::makeCluster(cores$cpus)")
		  cl = parallel::makeCluster(cores$cpus)
		}else{
		  print(paste("cores in use:",cores$cpus,sep=""))
		  print("cl = parallel::makeCluster(cores$cpus, type = 'FORK')")
		  cl = parallel::makeCluster(cores$cpus, type = "FORK")
		}

    #########################################################################
		#2016-11-26 begin
		#bbbb =	parallel::parLapply(cl, M, cooccur.gennetework.calculateCooccur, nrow, alpha, sequences$dt_idxtable, sequences$freqMatrix, sequences$matrix, sequences$constantList$biseqlevel, sequences$constantList$biseqidlevel)
		df_cooccurrence =	parallel::parSapply(cl, M, cooccur.gennetework.calculateCooccur, nrow, alpha, sequences$dt_idxtable, sequences$freqMatrix, sequences$matrix, sequences$constantList$biseqlevel, sequences$constantList$biseqidlevel)
		#print(df_cooccurrence)
		#########################################################################
		#2016-11-26 begin
		##bbbb <-	lapply(M, cooccur.gennetework.calculateCooccur, pb, progress, nrow, alpha, sequences$dt_idxtable, sequences$freqMatrix, sequences$matrix, sequences$constantList$biseqlevel, sequences$constantList$biseqidlevel, df_cooccurrence)
		parallel::stopCluster(cl)

		cooccur.printTimeCost("cooccur.gennetework.cooccurnetworks parallel::parLapply time cost",t,debug)
    #print(bbbb) #list

		t = as.character(Sys.time())
		#df_cooccurrence = Matrix::Matrix(nrow=nrow,ncol=0, sparse=TRUE, data=0)

		cooccur.printTimeCost("cooccur.gennetework.cooccurnetworks Matrix::cBind time cost",t,debug)
    #rm(bbbb)
    gc()
    #print(dim(df_cooccurrence))
    #print(debug)

		return(df_cooccurrence)
	}else if(sequences$memory=="sparssse" && parallel==TRUE){
	  #tttt = as.character(Sys.time())
	  #print(alpha)
	  #nrow = nrow(sequences$matrix)
	  #df_cooccurrence = sparseMatrix(nrow, nrow(sequences$dt_idxtable), x=0)
	  #df_cooccurrence = Matrix::Matrix(nrow=nrow,ncol=0,sparse=TRUE)
	  #if(!is.na(steps)){
	    #steps = steps+1
	    #cat(sprintf("%s. calculating networks ......", steps))
	    #message("")
	  #}

	  #print("foreach dopar ")
	  #pb <- txtProgressBar(style = 3)
	  #progress = seq(0,1, 1/(nrow(sequences$dt_idxtable)-1))

	  #for(n in 1:nrow(sequences$dt_idxtable)){
	  #library(doParallel)
	  #library(foreach)
	  #library(Matrix)
	  #gc()
	  #cores <- cooccur.detectCores()
	  #if(cores$win==TRUE){
	    #print("cl = parallel::makeCluster(cores$cpus)")
	    #cl = parallel::makeCluster(cores$cpus)
	  #}else{
	    #print("cl = parallel::makeCluster(cores$cpus, type = 'FORK')")
	    #cl = parallel::makeCluster(cores$cpus, "FORK")
	  #}
	  #registerDoParallel(cl);
	  #df_cooccurrence = foreach(n =  1:nrow(sequences$dt_idxtable), .packages='Matrix', .combine=cBind) %do% {
	  #df_cooccurrence = foreach(n =  1:nrow(sequences$dt_idxtable), .packages='Matrix', .combine=cBind) %dopar% {
	    #tttt = as.character(Sys.time())
	    #print(n)
	    #rowx = sequences$dt_idxtable[n,]
	    #i <- rowx[2]
	    #j <- rowx[3]
	    #i <- as.numeric(rowx[2])
	    #j <- as.numeric(rowx[3])
	    #x = round(sequences$freqMatrix[,i]/nrow,5)
	    #y = round(sequences$freqMatrix[,j]/nrow,5)
	    #xMy = x * y
	    #xy = sqrt(alpha * xMy)

	    #print(n)
	    #a = sequences$matrix[,i]
	    #b = sequences$matrix[,j]
	    #xx =   (round( a / b, 5) + a ) * 100000
	    #xx =    a * b  + a
	    #xx =   (round( sequences$matrix[,i] / sequences$matrix[,j],3) + sequences$matrix[,i]) * 1000
	    #xx = sequences$constantList$biseqlevel[match(xx, sequences$constantList$biseqidlevel)]

	    #print(xx)
	    #print(n)
	    #aa = table(xx)[xx]

	    #aa = as.vector(round(aa/nrow,5)	)
	    #bb = (which(aa>=xy))

	    #print(bb)
	    #print(n)
	    #m2 =  Matrix::Matrix(0,nrow=nrow,ncol=1,sparse=TRUE)
	    #if(length(bb)>=1){
	      #m2[bb,1]= 1
	    #}

	    #cooccur.printTimeCost('foreach test ',tttt,debug)


	    #return(m2)
	    #print(m2[,])

	    #df_cooccurrence <-  Matrix::cBind(df_cooccurrence, m2)
	    #setTxtProgressBar(pb, progress[n])
	  #}
	  #close(pb)
	  #stopCluster(cl);
	  #cooccur.printTimeCost('cooccur.gennetework.cooccurnetworks time cost',tttt,debug)
	  #print(df_cooccurrence)
	  #return(df_cooccurrence)
	}else if(sequences$memory=="sparse" && parallel==FALSE){
	  #print(alpha)
	  nrow = nrow(sequences$matrix)
	  #df_cooccurrence = sparseMatrix(nrow, nrow(sequences$dt_idxtable), x=0)
	  df_cooccurrence = Matrix::Matrix(nrow=nrow,ncol=0,sparse=TRUE)
	  if(!is.na(steps)){
	    steps = steps+1
	    cat(sprintf("%s. calculating networks ......", steps))
	    message("")
	  }

	  len = nrow(sequences$dt_idxtable)
	  colsperIter = 50
	  Iter = ceiling(len / colsperIter)
	  start <- 1
	  end <- 1

	  pb <- c()
	  if(debug){
  	  pb <- txtProgressBar(style = 3)
  	  if(Iter == 1){
  	    progress = seq(1, 1)
  	  }else{
  	    progress = seq(0, 1, 1/(Iter-1))
  	  }
    }

	  k = 1
	  for(k in 1:Iter){
	  #foreach::foreach(k =  1:Iter) %do% {
	    end <- k * colsperIter
	    if(end > len) end = len

	    i = sequences$dt_idxtable[,2][start:end]
	    j = sequences$dt_idxtable[,3][start:end]

	    x = round(sequences$freqMatrix[,i]/nrow,5)
	    y = round(sequences$freqMatrix[,j]/nrow,5)
	    xMy = x * y
	    xy = sqrt(alpha * xMy)


	    a = sequences$matrix[,i]
	    b = sequences$matrix[,j]
	    xx =   (round( a / b, 5) + a ) * 100000

	    #print(is.vector(xx))
      if(is.vector(xx)==FALSE){
  	    xx = apply(xx,2, function(x){sequences$constantList$biseqlevel[match(x, sequences$constantList$biseqidlevel)]})
  	    xx = (apply(xx, 2, function(x) {table(x)[x]}))
	    }else{
	      xx = sequences$constantList$biseqlevel[match(xx, sequences$constantList$biseqidlevel)]
	      xx = table(xx)[xx]
	      xx = as.matrix(xx)
	      #print(xx)
	    }
      #print(xx)
      #print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
	    xx = round(xx/nrow,5)

      #print(xy)
	    xx = ifelse(xx>=xy,1,0)
      #print(xx)

	    df_cooccurrence =  Matrix::cBind(df_cooccurrence, xx)
	    #print(df_cooccurrence)

      rm(xx)
      rm(xy)

	    start <- end + 1
	    end <- k * colsperIter

	    if(debug){
	      setTxtProgressBar(pb, progress[k])
	    }
	  }

    gc()
    if(debug){
	    close(pb)
    }
	  cooccur.printTimeCost('cooccur.gennetework.cooccurnetworks time cost',t,debug)
	  #print(df_cooccurrence)
	  #coocur.gennetwork.binom.test(df_cooccurrence)
	  return(df_cooccurrence)
	}
}

#2016-10-18
#' @importFrom stats binom.test
coocur.gennetwork.binom.test <- function(df_cooccurrence){
  #print(dim(df_cooccurrence))
  binomtest = apply(df_cooccurrence,2,function(x){
    if(sum(x)!=0){
      #xx = binom.test(sum(x>nrow(df_cooccurrence)),length(x),al="l");
      xx = binom.test(length(which(x==1)),length(x),p=round(length(which(x==1))/length(x),1),alternative="greater", conf.level = 0.99);
      #xx = binom.test(length(which(x==1)),length(x),p=round(length(which(x==1))/length(x),1), alternative="greater", conf.level = 0.99);

      #xx = ks.test(x,"pexp",0.01)
      #xx=wilcox.test(x)
      #print(xx)

      pvalue = round(xx$p.value,9)
      if(pvalue<0.01){
        pvalue="<0.01"
      }else if(pvalue==1){
        pvalue="1"
      }else{
        pvalue=round(pvalue,3)
      }
      #pvalue = round(xx$p.value,3)
      return(pvalue);
    }else{
      return(1)
    }
    #print(length(pvalue))
  })
  #print(which(binomtest<1))
  #print(length(binomtest))
  return(binomtest)
}

#2016-09-19
cooccur.gennetework.calculateCooccur  <- function(n, nrow, alpha,  dt_idxtable, freqMatrix, matrix, biseqlevel, biseqidlevel){
    #message("asdfasdfasd")
    #print(n)
    #browser()
    t = as.character(Sys.time())
    rowx = dt_idxtable[n,]
    i = rowx[2]
    j = rowx[3]
    #print(rowx)
    #i <- as.numeric(rowx[2])
    #j <- as.numeric(rowx[3])
    x = round(freqMatrix[,i]/nrow,5)
    y = round(freqMatrix[,j]/nrow,5)
    xy = x * y
    xy = sqrt(alpha * xy)


    a = matrix[,i]
    b = matrix[,j]
    xx =   (round( a / b, 5) + a ) * 100000
    #xx =    a * b  + a
    #xx =   (round( sequences$matrix[,i] / sequences$matrix[,j],3) + sequences$matrix[,i]) * 1000
    #2016-09-26 begin
    xx = biseqlevel[match(xx, biseqidlevel)]
    #xx = as.character(xx)

    #2016-09-26 end
    #print(xx)

    xx = table(xx)[xx]
    xx = as.vector(round(xx/nrow,5)	)
    bb = (which(xx>=xy))

    #print(bb)

    m2 =  Matrix::Matrix(nrow=nrow,ncol=1,sparse=TRUE,data=0)
    if(length(bb)>=1){
      #t = Sys.time()
      #df_cooccurrence[bb,n] = 1

      m2[bb,1]= 1


      #print(length(bb))
      #cooccur.printTimeCost('insert into df_cooccurrence',t,debug)
    }

    rm(a)
    rm(b)
    rm(xx)
    rm(bb)
    rm(xy)
    cooccur.printTimeCost('cooccur.gennetework.calculateCooccur ',t,TRUE)
    #print(m2[,])

    #2016-11-26 begin
    #return(m2)
    return(m2[,1])
}

cooccur.gennetework.calucation <- function(rowx, alpha=0.9,sequences){
	#print(rowx)



	nrow = nrow(sequences$matrix)


	#2016-11-16 begin
	#cooccurrenceList <- c()
	corrrry <- c(NA)
	length(corrrry) <- nrow
	#2016-11-16 end

	#i <- as.numeric(rowx["c"])

	#2016-11-20 begin
	#i = as.numeric(rowx[2])
	#j = as.numeric(rowx[3])

	#2016-11-16 begin
	##x = round(sequences$freqMatrix[,i]/nrow,5)
	##y = round(sequences$freqMatrix[,j]/nrow,5)
	##xMy = x * y

	#xMy = round(sequences$freqMatrix[,i]/nrow,5) * round(sequences$freqMatrix[,j]/nrow,5)
	xMy = round(sequences$freqMatrix[,as.numeric(rowx[2])]/nrow,5) * round(sequences$freqMatrix[,as.numeric(rowx[3])]/nrow,5)
	#2016-11-16 end

	#idx = which(sequences$dt_idxtable[,2]==i & sequences$dt_idxtable[,3]==j)
	idx = which(sequences$dt_idxtable[,2]==as.numeric(rowx[2]) & sequences$dt_idxtable[,3]==as.numeric(rowx[3]))
	#2016-11-20 end

	if(is.vector(sequences$bigramFreqList)){
		xy = round(sequences$bigramFreqList[idx]/nrow,5)
		#print(xy)
		#print(sqrt(alpha*x*y))
		corrrr = (unlist(xy))^2/xMy
	}else{
		xy = round(sequences$bigramFreqList[idx,]/nrow,5)
		#print(xy)
		#print(sqrt(alpha*x*y))
		corrrr = ((xy))^2/xMy
	}

	#print(corrrr)
	corrrry = ifelse(corrrr>=alpha, 1, 0)

	#2016-11-16 begin
	#cooccurrenceList= c(cooccurrenceList, corrrry)
	#print(cooccurrenceList)
	#return(cooccurrenceList)

	return(corrrry)
	#2016-11-16 end


}


cooccur.gennetework.outputNetWork.old <- function(i,sequences,df_cooccurrence){
	x = df_cooccurrence[i,]
	network <- list()
	network$Name = sequences$xnames[i]
	#cooccurrenceEdges <- c()
	#print(x)
	if(length(which(x==1))>0){

		sub_dfcooc = subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]%in%which(x==1))


		if(is.matrix(sub_dfcooc)){
			start = colnames(sequences$matrix)[sub_dfcooc[,2]]
			end = colnames(sequences$matrix)[sub_dfcooc[,3]]
		}else{
			start = colnames(sequences$matrix)[sub_dfcooc[2]]
			end = colnames(sequences$matrix)[sub_dfcooc[3]]
		}

		#cooccurrenceEdges = paste("(",start,",",end,")",sep="")
		network$Edges = paste(start,"-",end,sep="")

	}

	return(network)
}


cooccur.gennetework.outputNetWork <- function(i, sequences, df_cooccurrence, rawfile){
  x = df_cooccurrence[i,]
  network <- c()
  networkName = sequences$xnames[i]
  #cooccurrenceEdges <- c()
  #print(x)
  if(length(which(x==1))>0){

    sub_dfcooc = subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]%in%which(x==1))


    if(is.matrix(sub_dfcooc)){
      start = colnames(sequences$matrix)[sub_dfcooc[,2]]
      end = colnames(sequences$matrix)[sub_dfcooc[,3]]
    }else{
      start = colnames(sequences$matrix)[sub_dfcooc[2]]
      end = colnames(sequences$matrix)[sub_dfcooc[3]]
    }

    #cooccurrenceEdges = paste("(",start,",",end,")",sep="")
    network  = c(networkName, paste(start,"-",end,sep="", collapse = " "))
    write.table(paste(network,collapse =  " "), file=rawfile, append = TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
  }

  return(i)
}


cooccur.gennetework.outputNetWorkpvalue <- function(i, sequences, df_cooccurrence, networkpfile){
  x = df_cooccurrence[i,]
  network <- c()
  networkName = sequences$xnames[i]
  #cooccurrenceEdges <- c()
  #print(x)
  if(length(which(x>0))>0){

    sub_dfcooc = subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]%in%which(x>0))
    pvalues = round(x[which(x>0)],5)

    if(is.matrix(sub_dfcooc)){
      start = colnames(sequences$matrix)[sub_dfcooc[,2]]
      end = colnames(sequences$matrix)[sub_dfcooc[,3]]
    }else{
      start = colnames(sequences$matrix)[sub_dfcooc[2]]
      end = colnames(sequences$matrix)[sub_dfcooc[3]]
    }

    #cooccurrenceEdges = paste("(",start,",",end,")",sep="")
    network  = c(networkName, paste(start,"-",end,"-",pvalues,sep="", collapse = " "))
    write.table(paste(network,collapse =  " "), file=networkpfile, append = TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
  }

  return(i)
}

cooccur.gennetework.outputNewNetWorkFile <- function(i, sequences, df_cooccurrence, newnetworkfile){
  x = df_cooccurrence[i,]
  network <- c()
  networkName = sequences$xnames[i]
  #cooccurrenceEdges <- c()
  #print(x)
  if(length(which(x>0 & x<=0.05))>0){

    sub_dfcooc = subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]%in%which(x>0 & x<=0.05))
    pvalues = round(x[which(x>0 & x<=0.05)],5)

    if(is.matrix(sub_dfcooc)){
      start = colnames(sequences$matrix)[sub_dfcooc[,2]]
      end = colnames(sequences$matrix)[sub_dfcooc[,3]]
    }else{
      start = colnames(sequences$matrix)[sub_dfcooc[2]]
      end = colnames(sequences$matrix)[sub_dfcooc[3]]
    }

    #cooccurrenceEdges = paste("(",start,",",end,")",sep="")
    network  = c(networkName, paste(start,"-",end,"-",pvalues,sep="", collapse = " "))
    write.table(paste(network,collapse =  " "), file=newnetworkfile, append = TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
  }

  return(i)
}

cooccur.gennetework.outputNetWork.new <- function(x,dt_idxtable,matrix){
  #x = df_cooccurrence[i,]
  #network <- list()
  #network$Name = sequences$xnames[i]
  #cooccurrenceEdges <- c()
  #print(x)
  network = c()
  if(length(which(x==1))>0){

    sub_dfcooc = subset(dt_idxtable, dt_idxtable[,1]%in%which(x==1))


    if(is.matrix(sub_dfcooc)){
      start = colnames(matrix)[sub_dfcooc[,2]]
      end = colnames(matrix)[sub_dfcooc[,3]]
    }else{
      start = colnames(matrix)[sub_dfcooc[2]]
      end = colnames(matrix)[sub_dfcooc[3]]
    }

    #cooccurrenceEdges = paste("(",start,",",end,")",sep="")
    #network$Edges = paste(start,"-",end,sep="")
    network  = c(1, paste(start,"-",end,sep="", collapse = " "))

  }
  return(network)
}


cooccur.gennetework.outputNetWorkModule <- function(rowid,sequences,df_cooccurrence){

	x = df_cooccurrence[rowid,]

	moduleList <- list()
	moduleList$Name = sequences$xnames[rowid]

	module_str = ""
	#print(x)
	if(length(which(x==1))>0){
		#sub_dfcooc <- subset(sequences$dt_idxtable, b%in%which(x==TRUE))
		sub_dfcooc = subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]%in%which(x==1))
		#print(sub_dfcooc)


		#start <- colnames(sequences$matrix)[sub_dfcooc$c]
		#end <- colnames(sequences$matrix)[sub_dfcooc$d]

		if(is.matrix(sub_dfcooc)){
			start = colnames(sequences$matrix)[sub_dfcooc[,2]]
			end = colnames(sequences$matrix)[sub_dfcooc[,3]]
		}else{
			start = colnames(sequences$matrix)[sub_dfcooc[2]]
			end = colnames(sequences$matrix)[sub_dfcooc[3]]
		}
		#print(start)
		#print(end)

		module = c()
		for(i in 1:length(start)){
		#i = 1
		#foreach::foreach(i =  1:length(start)) %do% {
			s = start[i]
			e = end[i]
			#print(sprintf("index:%s,start:%s,end:%s",i,s,e))

			flag = FALSE
			if(length(module)>0){
				for(j in 1:length(module)){
				#j = 1
				#foreach::foreach(j =  1:length(module)) %do% {
					if(length(which(module[[j]]==s))>0 | length(which(module[[j]]==e))>0 ){
						flag = TRUE
						#print("found")
						if(length(which(module[[j]]==e))==0){
							module[[j]] = c(module[[j]], e)
						}
						if(length(which(module[[j]]==s))==0){
							module[[j]] = c(module[[j]], s)
						}
						break;
					}
					#Sys.sleep(5)
				}
			}
			if(flag == FALSE){
				#print("not found")
				submodule = c(s,e)
				module = c(module, list(submodule))
				#Sys.sleep(5)
			}
			#print(module)
		}

		if(length(module)>0){

			for(j in 1:length(module)){
			#j = 1
			#foreach::foreach(j =  1:length(module)) %do% {
				str = paste("(", paste(module[[j]],collapse = ","), ")",sep="")
				if(j==1){
					module_str = str
				}else{
					module_str = paste(module_str,str,sep=" ")
				}

			}
		}
		#print(module_str)

	}
	moduleList$Str = module_str
	return(moduleList)
}

cooccur.gennetework.outputNetWorkProperties <- function(rowid,sequences,df_cooccurrence){
	#require(igraph)
	if(!requireNamespace("igraph", quietly = TRUE)){
		stop("Package 'igraph' is required.")
	}
  #print(rowid)
	x = df_cooccurrence[rowid,]
	#print(x)
	properties = list()
	properties$Name = sequences$xnames[rowid]
	#print(sequences$xnames[rowid])
	if(length(which(x==1))>0){
		#print(length(which(x==1)))
		#print(ncol(sequences$original))
		################Connectivity######################
		connectivity = (length(which(x==1)) / sequences$original_ncol)
		if(length(connectivity)==0){
		  properties$Connectivity = 0
		}else{
		  properties$Connectivity = round(connectivity,6)
		}

		################################################
		###############prepare df#######################
		sub_dfcooc <- subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]%in%which(x==1))
		if(is.matrix(sub_dfcooc)){
			start = colnames(sequences$matrix)[sub_dfcooc[,2]]
			end = colnames(sequences$matrix)[sub_dfcooc[,3]]
		}else{
			start = colnames(sequences$matrix)[sub_dfcooc[2]]
			end = colnames(sequences$matrix)[sub_dfcooc[3]]
		}
		df = data.frame(A= start, B=end)
		########################################
		#################Diameter################
		df.g = igraph::graph.data.frame(d = df, directed = FALSE)
		properties$Diameter = igraph::diameter(df.g, directed = FALSE)
		########################################
		#################Radius################
		properties$Radius = igraph::radius(df.g)
		########################################
		################Entropy################
		#sigma_k = (length(which(x==1))) * 2
		#entropy = 0.0
		#print(unique(start))
		#edges = unique(c(start,end))
		#for(i in 1:length(edges)){
			#ss = edges[i]
			#len = length(which(start==ss)) + length(which(end==ss))
			#importance_i = len / sigma_k
			#print(len)
			#entropy = -importance_i * log(importance_i) + entropy;
		#}
		#properties$Entropy = round(entropy,6)
		########################################
		################ConnectionEffcient######
		#######clustering coefficient###########
		##2016-08-31
		properties$clusteringCoefficient = round(igraph::transitivity(df.g),6)
		########################################


	}
	#print(properties)
	return(properties)
}

#' @importFrom   foreach %do%
cooccur.gennetework.outputNetWorkcooccurrence.old <- function(sequences,df_cooccurrence,shuffeld=FALSE, debug){
	#require(Matrix)
	if(!requireNamespace("foreach", quietly = TRUE)){
		stop("Package 'foreach' is required.")
	}
	if(!requireNamespace("Matrix", quietly = TRUE)){
		stop("Package 'Matrix' is required.")
	}
	t = Sys.time()
	ncol = 0
	if(shuffeld==FALSE){
		ncol = 3
	}else{
		ncol = 1
	}
	cooccurrence = Matrix::Matrix(nrow=0,ncol=ncol,sparse=TRUE)

	cooccur.printTimeCost('Matrix::Matrix(nrow=0,ncol=ncol,sparse=TRUE) time cost',t,debug)

	#cooccurrence = Matrix(nrow=ncol(df_cooccurrence),ncol=ncol,sparse=TRUE)
	#for(colid in 1:ncol(df_cooccurrence)){
	colid =  1
	foreach::foreach(colid =  1:ncol(df_cooccurrence)) %do% {
		x = df_cooccurrence[,colid]
		#cooccurrence = list()
		cooccurr = c()
		#print(x)
		len = length(which(x==1))


		if(shuffeld==FALSE){
			sub_dfcooc = subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]==colid)
			#print(sub_dfcooc)
			#2016-11-20
			#start <- as.numeric(colnames(sequences$matrix)[sub_dfcooc[2]])
			#end <- as.numeric(colnames(sequences$matrix)[sub_dfcooc[3]])
			##cooccurrence$Site_i = start
			##cooccurrence$Site_j = end
			#cooccurr = c(start, end)
			cooccurr = c(as.numeric(colnames(sequences$matrix)[sub_dfcooc[2]]), as.numeric(colnames(sequences$matrix)[sub_dfcooc[3]]))
		}
		value = 0
		if(len>0){
			#print(len)
			#print(colid)
			value = round(len/length(x),3)
		}else{
			value = 0.000
		}
		cooccurr = c(cooccurr, value)
		cooccurrence = Matrix::rBind(cooccurrence,cooccurr)
		rm(cooccurr)
		gc()
		#cooccur.printTimeCost('foreach::foreach(colid =  1:ncol(df_cooccurrence)) time cost',t,debug)
		#cooccurrence[colid,] = cooccurr
	}
	cooccur.printTimeCost('cooccur.gennetework.outputNetWorkcooccurrence time cost',t,debug)
	return(cooccurrence)
}


#' @importFrom   foreach %do%
cooccur.gennetework.outputNetWorkcooccurrence <- function(sequences,df_cooccurrence,shuffeld=FALSE, parallel=FALSE, debug=FALSE){
  #require(Matrix)
  if(!requireNamespace("foreach", quietly = TRUE)){
    stop("Package 'foreach' is required.")
  }
  if(!requireNamespace("Matrix", quietly = TRUE)){
    stop("Package 'Matrix' is required.")
  }
  t = Sys.time()
  ncol = 0
  if(shuffeld==FALSE){
    ncol = 3
  }else{
    ncol = 1
  }

  #print(df_cooccurrence)
  t = Sys.time()
  if(parallel == TRUE){
    #cooccurrence = Matrix(nrow=ncol(df_cooccurrence),ncol=ncol,sparse=TRUE)
    cooccurrence = Matrix::Matrix(nrow=ncol(df_cooccurrence),ncol=ncol,sparse=TRUE,data=0)
    #for(colid in 1:ncol(df_cooccurrence)){
    colid =  1
    foreach::foreach(colid =  1:ncol(df_cooccurrence)) %do% {
      x = df_cooccurrence[,colid]
      #cooccurrence = list()
      cooccurr = c()
      #print(colid)
      len = length(which(x==1))


      if(shuffeld==FALSE){
        sub_dfcooc <- subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]==colid)
        #print(sub_dfcooc)
        #2016-11-20
        #start <- as.numeric(colnames(sequences$matrix)[sub_dfcooc[2]])
        #end <- as.numeric(colnames(sequences$matrix)[sub_dfcooc[3]])
        ##cooccurrence$Site_i = start
        ##cooccurrence$Site_j = end
        #cooccurr = c(start, end)
        cooccurr = c(as.numeric(colnames(sequences$matrix)[sub_dfcooc[2]]), as.numeric(colnames(sequences$matrix)[sub_dfcooc[3]]))
      }
      value = 0
      if(len>0){
        #print(len)
        #print(colid)
        value = round(len/length(x),3)
      }else{
        value = 0.000
      }
      cooccurr = c(cooccurr, value)
      #cooccurrence = Matrix::rBind(cooccurrence,cooccurr)
      #print(cooccurr)
      cooccurrence[colid,] = cooccurr

      if(colid %% 1000 == 0){
        cooccur.printTimeCost("cooccurrence[,colid] = cooccurr time cost",t,debug)
      }

      #rm(cooccurr)
      #gc()
      #cooccur.printTimeCost('foreach::foreach(colid =  1:ncol(df_cooccurrence)) time cost',t,debug)
      #cooccurrence[colid,] = cooccurr
    }
    cooccur.printTimeCost('cooccur.gennetework.outputNetWorkcooccurrence foreach::foreach  time cost',t,debug)
  }else if(parallel == FALSE){
    #cooccurrence = Matrix(nrow=ncol(df_cooccurrence),ncol=ncol,sparse=TRUE)
    cooccurrence = Matrix::Matrix(nrow=ncol(df_cooccurrence),ncol=ncol,sparse=TRUE,data=0)
    #for(colid in 1:ncol(df_cooccurrence)){
    colid =  1
    for(colid in  1:ncol(df_cooccurrence)){
    #foreach::foreach(colid =  1:ncol(df_cooccurrence)) %do% {
      x = df_cooccurrence[,colid]
      #cooccurrence = list()
      cooccurr = c()
      #print(colid)
      len = length(which(x==1))


      if(shuffeld==FALSE){
        sub_dfcooc <- subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]==colid)
        #print(sub_dfcooc)
        #2016-11-20
        #start <- as.numeric(colnames(sequences$matrix)[sub_dfcooc[2]])
        #end <- as.numeric(colnames(sequences$matrix)[sub_dfcooc[3]])
        ##cooccurrence$Site_i = start
        ##cooccurrence$Site_j = end
        #cooccurr = c(start, end)
        cooccurr = c(as.numeric(colnames(sequences$matrix)[sub_dfcooc[2]]), as.numeric(colnames(sequences$matrix)[sub_dfcooc[3]]))
      }
      value = 0
      if(len>0){
        #print(len)
        #print(colid)
        value = round(len/length(x),3)
      }else{
        value = 0.000
      }
      cooccurr = c(cooccurr, value)
      #cooccurrence = Matrix::rBind(cooccurrence,cooccurr)
      #print(cooccurr)
      cooccurrence[colid,] = cooccurr

      if(colid %% 1000 == 0){
        cooccur.printTimeCost("cooccurrence[,colid] = cooccurr time cost",t,debug)
      }

      #rm(cooccurr)
      #gc()
      #cooccur.printTimeCost('foreach::foreach(colid =  1:ncol(df_cooccurrence)) time cost',t,debug)
      #cooccurrence[colid,] = cooccurr
    }
    cooccur.printTimeCost('cooccur.gennetework.outputNetWorkcooccurrence foreach::foreach  time cost',t,debug)
  }else if(parallel == "TRUExx"){
    cores <- cooccur.detectCores()
    if(cores$win==TRUE){
      print("cl = parallel::makeCluster(cores$cpus)")
      cl = parallel::makeCluster(cores$cpus)
    }else{
      print("cl = parallel::makeCluster(cores$cpus, type = 'FORK')")
      cl = parallel::makeCluster(cores$cpus, "FORK")
    }

    #library(doParallel)
    #library(foreach)
    #library(Matrix)

    #registerDoParallel(cl);
    #cooccurrence = foreach(colid =  1:ncol(df_cooccurrence), .packages='Matrix', .combine=rbind) %do% {
    #cooccurrence = foreach(colid =  1:ncol(df_cooccurrence), .packages='Matrix', .combine=rbind) %dopar% {
      #foreach::foreach(colid =  1:ncol(df_cooccurrence)) %do% {
      x = df_cooccurrence[,colid]
      #cooccurrence = list()
      cooccurr = c()
      #print(colid)
      len = length(which(x==1))


      if(shuffeld==FALSE){
        sub_dfcooc = subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]==colid)
        #print(sub_dfcooc)
        #2016-11-20
        #start <- as.numeric(colnames(sequences$matrix)[sub_dfcooc[2]])
        #end <- as.numeric(colnames(sequences$matrix)[sub_dfcooc[3]])
        ##cooccurrence$Site_i = start
        ##cooccurrence$Site_j = end
        #cooccurr = c(start, end)
        cooccurr = c(as.numeric(colnames(sequences$matrix)[sub_dfcooc[2]]), as.numeric(colnames(sequences$matrix)[sub_dfcooc[3]]))
      }
      value = 0
      if(len>0){
        #print(len)
        #print(colid)
        value = round(len/length(x),3)
      }else{
        value = 0.000
      }
      cooccurr = c(cooccurr, value)
      #return(Matrix(cooccurr))
      return(cooccurr)

    #}
    #stopCluster(cl);
    #print(cooccurrence)
    cooccurrence <- Matrix::Matrix(cooccurrence)
    #return(cooccurrence)
    cooccur.printTimeCost('cooccur.gennetework.outputNetWorkcooccurrence foreach::foreach dopar time cost',t,debug)
  }

  cooccur.printTimeCost('cooccur.gennetework.outputNetWorkcooccurrence time cost',t,debug)

  #print(cooccurrence)
  return(cooccurrence)
}

cooccur.gennetework.genecooccurence <- function(x, dt_idxtable, matrix, shuffeld, debug){
  #x = df_cooccurrence[,colid]
  #cooccurrence = list()
  cooccurr = c()
  #print(x)
  len = length(which(x==1))


  if(shuffeld==FALSE){
    #sub_dfcooc <- subset(dt_idxtable, dt_idxtable[,1]==colid)
    #print(sub_dfcooc)
    #start <- as.numeric(colnames(matrix)[sub_dfcooc[2]])
    #end <- as.numeric(colnames(matrix)[sub_dfcooc[3]])
    #cooccurrence$Site_i = start
    #cooccurrence$Site_j = end
    cooccurr = c(0, 1)
  }
  value = 0
  if(len>0){
    #print(len)
    #print(colid)
    value = round(len/length(x),3)
  }else{
    value = 0.000
  }
  cooccurr = c(cooccurr, value)
  #print(cooccurr)
  #cooccurrence = Matrix::rBind(cooccurrence,cooccurr)
  #rm(cooccurr)
  #gc()
  return(cooccurr)
}



cooccur.gennetework.outputNetWorkcooccurrence.old <- function(colid,sequences,df_cooccurrence){
	x = df_cooccurrence[,colid]
	#cooccurrence = list()
	cooccurrence = c()
	#print(x)
	len = length(which(x==1))
	if(len>0){
		sub_dfcooc = subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]==colid)
		#print(sub_dfcooc)
		#start <- colnames(sequences$matrix)[sub_dfcooc[2]]
		#end <- colnames(sequences$matrix)[sub_dfcooc[3]]
		##cooccurrence$Site_i = start
		##cooccurrence$Site_j = end
		cooccurrence = c(cooccurrence, colnames(sequences$matrix)[sub_dfcooc[2]], colnames(sequences$matrix)[sub_dfcooc[3]])
		Cooccur = 0
		#if(len>0){
			Cooccur = round(len/length(x),3)
		#}else{
			#Cooccur = 0.000
		#}
		cooccurrence = as.numeric(c(cooccurrence, Cooccur))
	}
	return(cooccurrence)
}
###
