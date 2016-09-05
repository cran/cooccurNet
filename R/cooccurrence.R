
#' @title readseq
#'
#' @description Reading a sequence file in fasta format
#' @usage readseq(dataFile="",  dataType="DNA", debug=FALSE)
#' @examples
#' data = readseq(dataFile=getexample(dataType="protein"), dataType="protein")
#' @param dataFile  the file name of data
#' @param dataType  the type of data. It could be 'DNA' (default), 'protein', 'SNP' or 'other'
#' @param debug   FALSE by default; to indicate whether the debug message will be displayed or not
#' @return    a list object containing the original data matrix
#' @export
readseq <- function(dataFile="",  dataType="DNA", debug=FALSE){
  t = Sys.time()
  cat(("reading file......"))
  parallel=FALSE
  cpus=NA
  data <- cooccur.readfasta.read(dataFile=dataFile, parallel=parallel, cpus=cpus, debug=debug)
  data$dataType = dataType
  cat((sprintf("completed. (%s,%s)",nrow(data$matrix),ncol(data$matrix))))
  cooccur.printTimeCost('total co-occurrence.readsequence time cost',t,debug)
  return(data)
}



#' @title pprocess
#'
#' @description Filter the conservative columns (defined as the conservative score greater than the "conservativeFilter")
#' @usage pprocess(data=list(),  conservativeFilter=0.95, debug=FALSE)
#' @examples
#' data = readseq(dataFile=getexample(dataType="protein"),dataType="protein")
#' data_process = pprocess(data=data)
#' @param data   a list  return from function "readseq()"
#' @param conservativeFilter  0.95 by default. a number in the range of 0~1. The column with conservative score greater than it would be filtered in the later analyses;
#' @param debug   FALSE by default; to indicate whether the debug message will be displayed or not
#' @return  a list object containing the original data matrix and frequency matrix
#' @references Du, X., Wang, Z., Wu, A., Song, L., Cao, Y., Hang, H., & Jiang, T. (2008). Networks of genomic co-occurrence capture characteristics of human influenza A (H3N2) evolution. Genome research, 18(1), 178-187. doi:10.1101/gr.6969007
#' @export
pprocess <- function(data=list(),  conservativeFilter=0.95, debug=FALSE){
  parallel=FALSE
  cpus=NA
  t = Sys.time()
  cat("preprocess file......")
  dataType = data$dataType
  conservativeFilter = cooccur.conservativeFilter(conservativeFilter,dataType)
  
  threshold = conservativeFilter
  
  
  data <- cooccur.dataprepreprocess.preprocess(data=data, threshold=threshold, parallel=parallel, cpus=cpus, memory=NA, debug=debug )
  cat(sprintf("preprocess file......completed. (%s,%s)",nrow(data$matrix),ncol(data$matrix)))
  cooccur.printTimeCost('total co-occurrence.getco-occurrence time cost',t,debug)
  return(data)
}

#' @title gencooccur
#'
#' @description Construct the co-occurrence network
#' @examples
#' data = readseq(dataFile=getexample(dataType="protein"),dataType="protein")
#' data_process = pprocess(data=data)
#' #cooccurNetwork  = gencooccur(data=data_process)
#' @param data   a list return from function "pprocess()"
#' @param cooccurFilter  a number in the range of 0~1. It determines whether two columns are perfect co-occurrence;
#' @param networkFile   file name with full path for storing the co-occurrence network for each row;
#' @param module    FALSE by default. If it is set to be TRUE, the modules in each network of the networkFile would be calculated.
#' @param moduleFile   file name with full path  for storing the modules for co-occurrence network;
#' @param property   FALSE by default. If it is set to be TRUE, the properties for each network of the networkFile, including the network diameter, connectivity, ConnectionEffcient and so on, would be calculated.
#' @param propertyFile  character,  file name with full path  for storing the modules for co-occurrence network;
#' @param siteCo    FALSE by default. If it is set to be TRUE, the extent of co-occurrence between all pairs of columns would be calculated. It is defined as the ratio of rows with perfect co-occurrence.
#' @param siteCoFile  file name with full path  for storing the extent of co-occurrence between all pairs of columns, and the related p-values. The later are calculated by simulations as follows   firstly, all columns in the data are randomly permutated; then, the pairwise siteCos are calculated. This process would be repeated N times (the value depends on the parameter sampleTimes). For each pair of columns, the rank of the original siteCo in the N siteCos derived from simulations are considered as the p-value for the original siteCo.
#' @param sampleTimes   a integer of permutations in the simulation when calculating the p-values.
#' @param debug   FALSE by default; to indicate whether the debug message will be displayed or not
#' @return   a list and all the output file paths are attributed in it. \cr  The attribute "networkFile" stores the co-occurrence network for each row; \cr The attribute "moduleFile" is optional. When the module is set to be TRUE, it would be output. It stores the modules for co-occurrence network; \cr The attribute "propertyFile" is optional. When the property is set to be TRUE, it would be output. It stores the properties for co-occurrence network; \cr The attribute "siteCoFile" is optional. When the property is set to be TRUE, it would be output. It stores all the pairwise siteCos between columns.
#' @references Du, X., Wang, Z., Wu, A., Song, L., Cao, Y., Hang, H., & Jiang, T. (2008). Networks of genomic co-occurrence capture characteristics of human influenza A (H3N2) evolution. Genome research, 18(1), 178-187. doi:10.1101/gr.6969007
#' @export
gencooccur <- function(data=list(), cooccurFilter=NULL, networkFile='cooccurNetwork', module=FALSE,moduleFile='cooccurNetworkModule', property=FALSE, propertyFile='cooccurNetworkProperty', siteCo=FALSE, siteCoFile='siteCooccurr', sampleTimes=100, debug=FALSE){
  t = Sys.time()
  #sequences <- cooccur.gennetework.calculateNetWork(sequences, alpha=alpha, parallel=parallel, cpus=cpus, filterfile=filterfile, rawfile=rawfile, modulefile=modulefile, propertyfile=propertyfile, cooccurfile=cooccurfile, pvaluefile=pvaluefile, ptimes=ptimes,debug=debug)
  #cooccur.printTimeCost('total co-occurrence.getco-occurrence time cost',t,debug)
  #return(sequences)
  parallel=FALSE
  cpus=NA
  dataType = data$dataType
  
  cooccurFilter = cooccur.cooccurFilter(cooccurFilter,dataType)
  
  alpha = cooccurFilter
  ptimes = sampleTimes
  
  filterfile = NA
  cooccurfile=NA
  
  
  rawfile = networkFile
  
  
  if(module==TRUE){
    modulefile = moduleFile
  }else{
    modulefile=NA
  }
  
  if(property==TRUE){
    propertyfile = propertyFile
  }else{
    propertyfile=NA
  }
  
  if(siteCo==TRUE){
    pvaluefile = siteCoFile
  }else{
    pvaluefile=NA
  }
  
  
  
  data <- cooccur.gennetework.calculateNetWork(sequences=data, alpha=alpha, parallel=parallel, cpus=cpus, filterfile=filterfile, rawfile=rawfile, modulefile=modulefile, propertyfile=propertyfile, cooccurfile=cooccurfile, pvaluefile=pvaluefile, ptimes=ptimes,debug=debug)
  cooccur.printTimeCost('total co-occurrence.createco-occurrence time cost',t,debug)
  return(data)
  
}




#cooc <- function(filename, threshold=0.9, alpha=0.9, memory=NA, filterfile=NA, rawfile=NA, modulefile=NA, propertyfile=NA, cooccurfile=NA, pvaluefile=NA, ptimes=100, parallel=FALSE, cpus=NA,debug=FALSE){
#t = Sys.time()
#cat(("reading file......"))
#sequences <- cooccur.readfasta.read(filename, debug)
#cat((sprintf("completed. (%s,%s)",nrow(sequences$matrix),ncol(sequences$matrix))))
#message("")
#cat("preprocess file......")
#sequences <- cooccur.dataprepreprocess.preprocess(sequences,threshold, memory, debug)
#cat(sprintf("preprocess file......completed. (%s,%s)",nrow(sequences$matrix),ncol(sequences$matrix)))
#message("")
#sequences <- cooccur.gennetework.calculateNetWork(sequences, alpha=alpha, parallel=parallel, cpus=cpus, filterfile=filterfile, rawfile=rawfile, modulefile=modulefile, propertyfile=propertyfile, cooccurfile=cooccurfile, pvaluefile=pvaluefile, ptimes=ptimes,debug=debug)
#cooccur.printTimeCost('total co-occurrence.createco-occurrence time cost',t,debug)
#return(sequences)
#}





#' @title coocnet
#'
#' @description Read and preprocess data, and construct the co-occurrence network in one step.
#' @examples
#' cooccurNetwork  = coocnet(dataFile=getexample(dataType="protein"), dataType="protein")
#' @param dataFile  file name with full path of DNA, protein, SNP data or other kinds of data
#' @param dataType  the type of data. It could be 'DNA' (default), 'protein', 'SNP' or 'other'
#' @param conservativeFilter  0.95 by default. a number in the range of 0~1. The column with conservative score greater than it would be filtered in the later analyses;
#' @param cooccurFilter  a number in the range of 0~1. It determines whether two columns are perfect co-occurrence;
#' @param networkFile   file name with full path for storing the co-occurrence network for each row;
#' @param module    FALSE by default. If it is set to be TRUE, the modules in each network of the networkFile would be calculated.
#' @param moduleFile   file name with full path  for storing the modules for co-occurrence network;
#' @param property   FALSE by default. If it is set to be TRUE, the properties for each network of the networkFile, including the network diameter, connectivity, ConnectionEffcient and so on, would be calculated.
#' @param propertyFile  character,  file name with full path  for storing the modules for co-occurrence network;
#' @param siteCo    FALSE by default. If it is set to be TRUE, the extent of co-occurrence between all pairs of columns would be calculated. It is defined as the ratio of rows with perfect co-occurrence.
#' @param siteCoFile  file name with full path  for storing the extent of co-occurrence between all pairs of columns, and the related p-values. The later are calculated by simulations as follows   firstly, all columns in the data are randomly permutated; then, the pairwise siteCos are calculated. This process would be repeated N times (the value depends on the parameter sampleTimes). For each pair of columns, the rank of the original siteCo in the N siteCos derived from simulations are considered as the p-value for the original siteCo.
#' @param sampleTimes   a integer of permutations in the simulation when calculating the p-values.
#' @param debug   FALSE by default; to indicate whether the debug message will be displayed or not
#' @return   a list and all the output file paths are attributed in it. \cr  The attribute "networkFile" stores the co-occurrence network for each row; \cr The attribute "moduleFile" is optional. When the module is set to be TRUE, it would be output. It stores the modules for co-occurrence network; \cr The attribute "propertyFile" is optional. When the property is set to be TRUE, it would be output. It stores the properties for co-occurrence network; \cr The attribute "siteCoFile" is optional. When the property is set to be TRUE, it would be output. It stores all the pairwise siteCos between columns.
#' @references Du, X., Wang, Z., Wu, A., Song, L., Cao, Y., Hang, H., & Jiang, T. (2008). Networks of genomic co-occurrence capture characteristics of human influenza A (H3N2) evolution. Genome research, 18(1), 178-187. doi:10.1101/gr.6969007
#' @export
coocnet <-  function(dataFile="", dataType="DNA", conservativeFilter=0.95, cooccurFilter=NULL, networkFile='cooccurNetwork', module=FALSE,moduleFile='cooccurNetworkModule', property=FALSE, propertyFile='cooccurNetworkProperty', siteCo=FALSE, siteCoFile='siteCooccurr', sampleTimes=100,  debug=FALSE){
  parallel=FALSE
  cpus=NA
  t = Sys.time()
  cat(("reading file......"))
  sequences <- cooccur.readfasta.read(dataFile=dataFile, debug=debug)
  cat((sprintf("completed. (%s,%s)",nrow(sequences$matrix),ncol(sequences$matrix))))
  message("")
  cat("preprocess file......")
  
  conservativeFilter = cooccur.conservativeFilter(conservativeFilter,dataType)
  
  cooccurFilter = cooccur.cooccurFilter(cooccurFilter,dataType)
  
  threshold = conservativeFilter
  memory=NA
  sequences <- cooccur.dataprepreprocess.preprocess(data=sequences,threshold=threshold, memory=memory, debug=debug)
  cat(sprintf("preprocess file......completed. (%s,%s)",nrow(sequences$matrix),ncol(sequences$matrix)))
  message("")
  
  alpha = cooccurFilter
  ptimes = sampleTimes
  
  filterfile = NA
  cooccurfile=NA
  
  
  rawfile = networkFile
  
  
  if(module==TRUE){
    modulefile = moduleFile
  }else{
    modulefile=NA
  }
  
  if(property==TRUE){
    propertyfile = propertyFile
  }else{
    propertyfile=NA
  }
  
  if(siteCo==TRUE){
    pvaluefile = siteCoFile
  }else{
    pvaluefile=NA
  }
  
  sequences <- cooccur.gennetework.calculateNetWork(sequences=sequences, alpha=alpha, parallel=parallel, cpus=cpus, filterfile=filterfile, rawfile=rawfile, modulefile=modulefile, propertyfile=propertyfile, cooccurfile=cooccurfile, pvaluefile=pvaluefile, ptimes=ptimes,debug=debug)
  cooccur.printTimeCost('total co-occurrence.createco-occurrence time cost',t,debug)
  return(sequences)
  
  
}




#' @title siteco
#'
#' @description Read and preprocess data, and calculate the pairwise site co-occurrence in one step
#' @examples
#' #pairwiseCooccur = siteco(dataFile=getexample(dataType="protein"), dataType="protein")
#' @param dataFile  file name with full path of DNA, protein, SNP data or other kinds of data
#' @param dataType  the type of data. It could be 'DNA' (default), 'protein', 'SNP' or 'other'
#' @param conservativeFilter 0.95 by default. a number in the range of 0~1. The column with conservative score greater than it would be filtered in the later analyses;
#' @param cooccurFilter  a number in the range of 0~1. It determines whether two columns are perfect co-occurrence;
#' @param siteCoFile  file name with full path  for storing the extent of co-occurrence between all pairs of columns, and the related p-values. The later are calculated by simulations as follows   firstly, all columns in the data are randomly permutated; then, the pairwise siteCos are calculated. This process would be repeated N times (the value depends on the parameter sampleTimes). For each pair of columns, the rank of the original siteCo in the N siteCos derived from simulations are considered as the p-value for the original siteCo.
#' @param sampleTimes   a integer of permutations in the simulation when calculating the p-values.
#' @param debug   FALSE by default; to indicate whether the debug message will be displayed or not
#' @return   a list and the output file path of 'sitecoFile' is attributed in it. The file stores all the pairwise siteCos between columns.
#' @references Du, X., Wang, Z., Wu, A., Song, L., Cao, Y., Hang, H., & Jiang, T. (2008). Networks of genomic co-occurrence capture characteristics of human influenza A (H3N2) evolution. Genome research, 18(1), 178-187. doi:10.1101/gr.6969007
#' @export
siteco <-  function(dataFile="", dataType="DNA", conservativeFilter=0.95, cooccurFilter=NULL, siteCoFile='siteCooccurr', sampleTimes=100, debug=FALSE){
  parallel=FALSE
  cpus=NA
  t = Sys.time()
  cat(("reading file......"))
  sequences <- cooccur.readfasta.read(dataFile=dataFile, debug=debug)
  cat((sprintf("completed. (%s,%s)",nrow(sequences$matrix),ncol(sequences$matrix))))
  message("")
  cat("preprocess file......")
  
  conservativeFilter = cooccur.conservativeFilter(conservativeFilter,dataType)
  
  cooccurFilter = cooccur.cooccurFilter(cooccurFilter,dataType)
  
  threshold = conservativeFilter
  memory=NA
  sequences <- cooccur.dataprepreprocess.preprocess(data=sequences,threshold=threshold, memory=memory, debug=debug)
  cat(sprintf("preprocess file......completed. (%s,%s)",nrow(sequences$matrix),ncol(sequences$matrix)))
  message("")
  
  alpha = cooccurFilter
  ptimes = sampleTimes
  
  filterfile = NA
  cooccurfile=NA
  rawfile=NA
  modulefile=NA
  propertyfile=NA
  pvaluefile = siteCoFile
  
  
  sequences <- cooccur.gennetework.calculateNetWork(sequences=sequences, alpha=alpha, parallel=parallel, cpus=cpus, filterfile=filterfile, rawfile=rawfile, modulefile=modulefile, propertyfile=propertyfile, cooccurfile=cooccurfile, pvaluefile=pvaluefile, ptimes=ptimes,debug=debug)
  cooccur.printTimeCost('total co-occurrence.createco-occurrence time cost',t,debug)
  return(sequences)
  
  
}

#' @title toigraph
#'
#' @description transform a network file to the igraph.data.frame type by specifying a network file name and a network name
#' @usage toigraph(networkFile="", networkNames=c())
#' @examples
#' #cooccurNetwork  = cooc(dataFile=getexample(dataType="protein"),dataType="protein")
#' #network_igraph = toigraph(networkFile=cooccurNetwork$networkFile, networkName=c("EPI823725")
#' @importFrom   stats na.omit
#' @importFrom   utils read.csv
#' @param networkFile the network file name for the co-occurrence network. It is generated from cooc() or gencooccur().
#' @param networkNames a vector of network names
#' @return   a list of igraph::graph.data.frame object ordered by the input network names.
#' @export
toigraph <-  function(networkFile="", networkNames=c()){
  
  if(!requireNamespace("igraph", quietly = TRUE)){
    stop("Package 'igraph' is required.")
  }
  
  dataframe = read.csv(networkFile,sep="\t",header=FALSE,stringsAsFactors = FALSE)
  if( is.null(dataframe)){
    stop("the file is empty")
  }
  if(nrow(dataframe)==0){
    stop("the file is empty")
  }
  
  if(length(networkNames)==0){
    networkNames = dataframe[,1]
  }
  
  graphList = list()
  
  for(i in 1:length(networkNames)){
    networkName = networkNames[i]
    
    df = dataframe[which(dataframe[,1]==networkName),]
    if(nrow(df)==0){
      #stop("the file does not contain the specified network.")
      graphList[[i]] = NA
      next
    }
    graphdf = matrix(ncol=3)
    for(rowid in  1:nrow(df)){
      rows = df[rowid,]
      NAME = as.character(rows[1])
      for (colid in  2:ncol(df)) {
        
        #newRow = c()
        cols = as.character(rows[colid])
        #print(NAME)
        #print(cols)
        if(cols!=""){
          value = strsplit(cols,split="-")
          A = as.numeric(value[[1]][1])
          B = as.numeric(value[[1]][2])
          #newRow =  c(NAME, A, B)
          #print(c(NAME, A, B))
          #graphdf =  rbind(c(NAME, A, B),graphdf)
          graphdf =  rbind(graphdf, c(NAME, A, B))
        }else{
          break
        }
      }
    }
    graphdf = na.omit(graphdf)
    
    graphdf = graphdf[which(graphdf[,1]==networkName),2:3]
    graphdf = igraph::graph.data.frame(d = graphdf, directed = FALSE)
    
    
    graphList[[i]] = graphdf
    
  }
  
  return(graphList)
}





#'getexample
#' @title getexample
#' @description Get the example data
#' @note  Both the dataset "DNA" and "protein" are sampled from the Hemagglutinin sequences of human influenza H3N2 viruses, while the dataset "SNP" are simulated.
#' @usage getexample(dataType)
#' @examples
#' dataFile = getexample(dataType='protein')
#' @param dataType  the type of data. It could be 'DNA' (default), 'protein', 'SNP' or 'other'
#' @return  file path of the example data
#' @export
#'
getexample <- function(dataType){
  path = ""
  if(dataType=="DNA"){
    path = system.file("extdata", "DNA", package = "cooccurNet")
  }else if(dataType=="protein"){
    path = system.file("extdata", "protein", package = "cooccurNet")
  }else if(dataType=="SNP"){
    path = system.file("extdata", "SNP", package = "cooccurNet")
  }else{
    stop("the param dataType should in the range of [DNA, protein,SNP]")
  }
  
  return(path)
}


#'getexample_forRCOS
#' @title getexample_forRCOS
#' @description Get the test data for testing RCOS method
#' @note  File paths for all the files available for testing RCOS. Currently, there is only one file "HA1protein_humanH3N2", the sequences within which are derived from the database of Influenza Virus Resource.
#' @usage getexample_forRCOS()
#' @examples
#' dataFile = getexample_forRCOS()
#' @return  File paths for all the files available for testing RCOS.
#' @export
#'
getexample_forRCOS <- function(){
  path = ""
  path = system.file("extdata", "RCOS", package = "cooccurNet")
  fileName<-paste(path, dir(path) ,sep="/")
  #print(path)
  return(fileName)
}
