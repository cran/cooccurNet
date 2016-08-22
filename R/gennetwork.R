cooccur.gennetework.calculateNetWork <- function(sequences=list(),alpha=0.9,parallel=FALSE,cpus=cooccur.detectCores(),filterfile=NA,rawfile=NA,modulefile=NA, propertyfile=NA, cooccurfile=NA, pvaluefile=NA, ptimes=100,debug=FALSE){

	steps = 0

	cooccurrenceList = c()
	cooccurrenceEdges = c()

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
		header$h = "name	Connectivity	Diameter	Radius	ConnectionEffcient"
		cooccur.writetable(header$h,propertyfile)
	}


	if(!is.na(pvaluefile)){
		sequences$siteCoFile = pvaluefile
		if(file.exists(pvaluefile)){
			file.remove(pvaluefile)
		}
		header = list()
		header$h = "Site_i	Site_j	Co-occur	p-value"
		cooccur.writetable(header$h,pvaluefile)
	}


	message("")

	df_cooccurrence = cooccur.gennetework.cooccurnetworks(sequences,alpha,steps)
	if(sequences$memory=="sparse") {
		steps = steps + 1
	}
	colnames(df_cooccurrence) =seq(1:ncol(df_cooccurrence))
	rownames(df_cooccurrence) = seq(1:nrow(sequences$matrix))


	if(!is.na(rawfile)){
		steps = steps+1
		cat(sprintf("%s. creating network file (%s) ......", steps,rawfile))
		t = Sys.time()
		cooccurrenceEdgesBM <- data.frame()
		cooccurrenceEdgesBM <- lapply(seq_len(nrow(df_cooccurrence)), cooccur.gennetework.outputNetWork, sequences, df_cooccurrence)
		cooccur.printTimeCost('create cooccurrence Edges time cost',t,debug)
		t = Sys.time()
		cooccur.writelist(cooccurrenceEdgesBM,rawfile)
		cooccur.printTimeCost('write networks networkFile time cost',t, debug)
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
		#t = Sys.time()
		#cooccurrence <- lapply(seq_len(ncol(df_cooccurrence)), cooccur.gennetework.outputNetWorkcooccurrence, sequences, df_cooccurrence)
		cooccurrence <- cooccur.gennetework.outputNetWorkcooccurrence(sequences, df_cooccurrence)
		#print(cooccurrence)
		#cooccur.printTimeCost('create cooccurrence cooccurrence time cost',t,debug)
		t = Sys.time()
		#cooccur.writelist(cooccurrence,cooccurfile)
		cooccur.writeMatrix(cooccurrence,"cooccurfile")
		cooccur.printTimeCost('write networks siteCoFile time cost',t,debug)
		cat("completed")
	}
	message("")
	if(!is.na(pvaluefile)){
		steps = steps+1
		cat(sprintf("%s. creating siteCoFile (%s), sampleTimes: %s ......",steps, pvaluefile,ptimes))


		#if(nrow(cooccurrence)==0){
		#if(length(cooccurrence)==0){
			#cooccurrence <- lapply(seq_len(ncol(df_cooccurrence)), cooccur.gennetework.outputNetWorkcooccurrence, sequences, df_cooccurrence)
			#print("-----------------------------------------------------------------------------------------")
			cooccurrence <- cooccur.gennetework.outputNetWorkcooccurrence(sequences, df_cooccurrence)
			#print(cooccurrence)
		#}#else{
			#cooccurrence = do.call("rbind", cooccurrence)
		#}

		colnames(cooccurrence) = c("Site_i","Site_j","Cooccur")
		t = Sys.time()
		pvalues = cooccur.networkpvalue.calculateNetWorkPvalue(cooccurrence, sequences, ptimes, alpha)
		cooccur.printTimeCost('create cooccurrence p-value time cost',t,debug)
		t = Sys.time()
		cooccur.writetable(pvalues,pvaluefile)
		cooccur.printTimeCost('write networks siteCoFile time cost',t,debug)
		cat("completed")
	}

	return(sequences)

}


#' @importFrom   foreach %do%
cooccur.gennetework.cooccurnetworks <- function(sequences, alpha=0.9, steps){
	if(!requireNamespace("foreach", quietly = TRUE)){
		stop("Package 'foreach' is required.")
	}
	t = Sys.time()
	if(sequences$memory=="memory"){

		df_cooccurrence <- data.frame()
		cooccurrenceList = c()
		for(i in 1:nrow(sequences$dt_idxtable)){
			rowx = sequences$dt_idxtable[i,]
			cooccurrence = cooccur.gennetework.calucation(rowx,alpha,sequences)
			cooccurrenceList = c(cooccurrenceList,list(cooccurrence))
		}
		df_cooccurrence <- t(do.call(rbind, cooccurrenceList))
		#cooccur.printTimeCost('cooccur.gennetework.cooccurnetworks time cost',t,debug)
		return(df_cooccurrence)
	}else if(sequences$memory=="sparse"){
		nrow = nrow(sequences$matrix)
		#df_cooccurrence = sparseMatrix(nrow, nrow(sequences$dt_idxtable), x=0)
		df_cooccurrence = Matrix::Matrix(nrow=nrow,ncol=0,sparse=TRUE)
		if(!is.na(steps)){
			steps = steps+1
			cat(sprintf("%s. calculating networks ......", steps))
			message("")
		}

		pb <- txtProgressBar(style = 3)
		progress = seq(0,1, 1/(nrow(sequences$dt_idxtable)-1))

		#for(n in 1:nrow(sequences$dt_idxtable)){
		n = 1
		foreach::foreach(n =  1:nrow(sequences$dt_idxtable)) %do% {
			#print(n)
			rowx = sequences$dt_idxtable[n,]
			i <- rowx[2]
			j <- rowx[3]
			#i <- as.numeric(rowx[2])
			#j <- as.numeric(rowx[3])
			x = round(sequences$freqMatrix[,i]/nrow,5)
			y = round(sequences$freqMatrix[,j]/nrow,5)
			xMy = x * y
			xy = sqrt(alpha * xMy)


			a = sequences$matrix[,i]
			b = sequences$matrix[,j]
			xx =   (round( a / b, 5) + a ) * 100000
			#xx =    a * b  + a 
			#xx =   (round( sequences$matrix[,i] / sequences$matrix[,j],3) + sequences$matrix[,i]) * 1000
			xx = sequences$constantList$biseqlevel[match(xx, sequences$constantList$biseqidlevel)]


			aa = table(xx)[xx]
			aa = as.vector(round(aa/nrow,5)	)
			bb = (which(aa>=xy))

			m2 =  Matrix::Matrix(0,nrow=nrow,ncol=1,sparse=TRUE)
			if(length(bb)>1){
				#t = Sys.time()
				#df_cooccurrence[bb,n] = 1

				m2[bb,1]= 1


				#print(length(bb))
				#cooccur.printTimeCost('insert into df_cooccurrence',t,debug)
			}
			df_cooccurrence <-  Matrix::cBind(df_cooccurrence, m2)
			setTxtProgressBar(pb, progress[n])
		}
		close(pb)
		#cooccur.printTimeCost('cooccur.gennetework.cooccurnetworks time cost',t,debug)
		return(df_cooccurrence)
	}
}


cooccur.gennetework.calucation <- function(rowx, alpha=0.9,sequences){
	#print(rowx)
	cooccurrenceList <- c()
	nrow = nrow(sequences$matrix)
	#i <- as.numeric(rowx["c"])
	i <- as.numeric(rowx[2])

	j <- as.numeric(rowx[3])

	x = round(sequences$freqMatrix[,i]/nrow,5)

	y = round(sequences$freqMatrix[,j]/nrow,5)

	xMy = x * y

	idx <- which(sequences$dt_idxtable[,2]==i & sequences$dt_idxtable[,3]==j)


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
	corrrry = ifelse(corrrr>=alpha, 1, 0)
	cooccurrenceList= c(cooccurrenceList, corrrry)
	return(cooccurrenceList)
}


cooccur.gennetework.outputNetWork <- function(i,sequences,df_cooccurrence){
	x = df_cooccurrence[i,]
	network <- list()
	network$Name = sequences$xnames[i]
	#cooccurrenceEdges <- c()
	#print(x)
	if(length(which(x==1))>0){

		sub_dfcooc <- subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]%in%which(x==1))


		if(is.matrix(sub_dfcooc)){
			start <- colnames(sequences$matrix)[sub_dfcooc[,2]]
			end <- colnames(sequences$matrix)[sub_dfcooc[,3]]
		}else{
			start <- colnames(sequences$matrix)[sub_dfcooc[2]]
			end <- colnames(sequences$matrix)[sub_dfcooc[3]]
		}

		#cooccurrenceEdges = paste("(",start,",",end,")",sep="")
		network$Edges = paste(start,"-",end,sep="")

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
		sub_dfcooc <- subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]%in%which(x==1))
		#print(sub_dfcooc)


		#start <- colnames(sequences$matrix)[sub_dfcooc$c]
		#end <- colnames(sequences$matrix)[sub_dfcooc$d]

		if(is.matrix(sub_dfcooc)){
			start <- colnames(sequences$matrix)[sub_dfcooc[,2]]
			end <- colnames(sequences$matrix)[sub_dfcooc[,3]]
		}else{
			start <- colnames(sequences$matrix)[sub_dfcooc[2]]
			end <- colnames(sequences$matrix)[sub_dfcooc[3]]
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
	properties <- list()
	properties$Name = sequences$xnames[rowid]
	#print(sequences$xnames[rowid])
	if(length(which(x==1))>0){
		#print(length(which(x==1)))
		#print(ncol(sequences$original))
		################Connectivity######################
		connectivity = (length(which(x==1)) / sequences$original_ncol)
		properties$Connectivity = round(connectivity,6)
		################################################
		###############prepare df#######################
		sub_dfcooc <- subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]%in%which(x==1))
		if(is.matrix(sub_dfcooc)){
			start <- colnames(sequences$matrix)[sub_dfcooc[,2]]
			end <- colnames(sequences$matrix)[sub_dfcooc[,3]]
		}else{
			start <- colnames(sequences$matrix)[sub_dfcooc[2]]
			end <- colnames(sequences$matrix)[sub_dfcooc[3]]
		}
		df <- data.frame(A= start, B=end)
		########################################
		#################Diameter################
		df.g <- igraph::graph.data.frame(d = df, directed = FALSE)
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
		properties$ConnectionEffcient = round(igraph::transitivity(df.g),6)
		########################################


	}
	return(properties)
}

#' @importFrom   foreach %do%
cooccur.gennetework.outputNetWorkcooccurrence <- function(sequences,df_cooccurrence,shuffeld=FALSE){
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
			sub_dfcooc <- subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]==colid)
			#print(sub_dfcooc)
			start <- as.numeric(colnames(sequences$matrix)[sub_dfcooc[2]])
			end <- as.numeric(colnames(sequences$matrix)[sub_dfcooc[3]])
			#cooccurrence$Site_i = start
			#cooccurrence$Site_j = end
			cooccurr = c(start, end)
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
		#cooccurrence[colid,] = cooccurr
	}
	#cooccur.printTimeCost('calculate network cooccurrence time cost',t,debug)
	return(cooccurrence)
}


cooccur.gennetework.outputNetWorkcooccurrence.old <- function(colid,sequences,df_cooccurrence){
	x = df_cooccurrence[,colid]
	#cooccurrence = list()
	cooccurrence = c()
	#print(x)
	len = length(which(x==1))
	if(len>0){
		sub_dfcooc <- subset(sequences$dt_idxtable, sequences$dt_idxtable[,1]==colid)
		#print(sub_dfcooc)
		start <- colnames(sequences$matrix)[sub_dfcooc[2]]
		end <- colnames(sequences$matrix)[sub_dfcooc[3]]
		#cooccurrence$Site_i = start
		#cooccurrence$Site_j = end
		cooccurrence = c(cooccurrence, start, end)
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
