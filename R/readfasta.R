


#' @importFrom   stats rnorm
cooccur.readfasta.read <- function(dataFile="", parallel=FALSE, cpus=cooccur.detectCores(), debug=FALSE){
  #require(seqinr)
  #require(data.table)
  #require(snowfall)

	if(!requireNamespace("seqinr", quietly = TRUE)){
		stop("Package 'seqinr' is required.")
	}

	if(!requireNamespace("data.table", quietly = TRUE)){
		stop("Package 'data.table' is required.")
	}

	#if(!requireNamespace("snowfall", quietly = TRUE)){
		#stop("Package 'snowfall' is required.")
	#}


	seqlevel = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","-","0","1","2","3","4","5","6","7","8","9")
	seqidlevel = c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,73, 79, 83, 89, 97, 101, 103, 107, 109, 113,127,131,137,139,149,151,157,163,173,179,181,191,193  )
	cooccur.preprocess.object = list(matrix=NA, freqMatrix=NA, dt_idxtable=NA, bigramFreqList=NA, memory=NA, conn=NA, seqlevel=NA, seqidlevel=NA)

	cooccur.preprocess.object$seqlevel = seqlevel
	cooccur.preprocess.object$seqidlevel = seqidlevel

	sequences = c()

	t = Sys.time()
	sequences <- data.table::fread(dataFile, header = FALSE, na.strings = "NA")
	sequences <-  sequences$V1

	cooccur.printTimeCost('reading sequence from file time cost',t,debug)
	#
	# Get the line numbers where sequences names are:
	#
	ind <- which(substr(sequences, 1L, 1L) == ">")

	xnames = sequences[ind]
	xnames = substr(xnames, 2, nchar(xnames))
	cooccur.preprocess.object$xnames = xnames

	nseq <- length(ind)
	if(nseq == 0){
		stop("no line starting with a > character found")
	}
	#
	# Localize sequence data:
	#
	start <- ind + 1
	end <- ind - 1
	end <- c(end[-1], length(sequences))
	#
	# Read sequences:
	#
	t = Sys.time()

	sequences <- lapply(seq_len(nseq), function(i) paste(sequences[start[i]:end[i]], collapse = ""))

	sequences_original <- sequences
	sequences <- lapply(sequences, function(x){
		x = seqinr::s2c(x)
		x = match(x, seqlevel)
		x = seqidlevel[x]
	})

	sequences_original <- lapply(sequences_original, function(x){
		x = seqinr::s2c(x)
	})

	sequences <- do.call(rbind, sequences)
	sequences_original <- do.call(rbind, sequences_original)

	colnames(sequences) <- seq(1:ncol(sequences))


	cooccur.preprocess.object$original = sequences_original
	cooccur.preprocess.object$original_ncol = ncol(sequences_original)

	cooccur.preprocess.object$matrix = sequences
	storage.mode(cooccur.preprocess.object$matrix) <- "integer"

	constantList <- cooccur.constant()
	cooccur.preprocess.object$constantList = constantList

	dt_idxtable_filename = paste("dt_idxtable_",rnorm(1)*100000000,sep="")
	#print(dt_idxtable_filename)
	cooccur.preprocess.object$dt_idxtable_filename = dt_idxtable_filename
	cooccur.printTimeCost('split string to character time cost',t, debug)
	return(cooccur.preprocess.object)
}




###
