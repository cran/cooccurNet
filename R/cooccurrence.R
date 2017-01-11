
#' @title readseq
#'
#' @description Reading a sequence file in fasta format
#' @usage readseq(dataFile="", dataType="protein", debug=FALSE)
#' @examples
#' data = readseq(dataFile=getexample(dataType="protein"), dataType="protein")
#' @param dataFile character, a FASTA data file name with full path.
#' @param dataType character, 'protein' by default, the type of data will be processed. It could be 'DNA', 'RNA', 'protein', 'SNP' or 'other'.
#' @param debug  logic, FALSE by default, to indicate whether the debug message will be displayed or not.
#' @return  list,contains the original data matrix and some other information.
#' @export
readseq <- function(dataFile="", dataType="protein", debug=FALSE){
  t = Sys.time()
  cat(("Reading file......"))
  #parallel=FALSE
  #cpus=NA
  data <- cooccur.readfasta.read(dataFile=dataFile, debug=debug)
  data$dataType = dataType
  cat((sprintf("completed. (%s,%s)",nrow(data$matrix),ncol(data$matrix))))
  cooccur.printTimeCost('total co-occurrence.readsequence time cost',t,debug)
  message("")
  return(data)
}



#' @title pprocess
#'
#' @description Filter the conservative columns (defined as the conservative score greater than the "conservativeFilter")
#' @usage pprocess(data=list(),  conservativeFilter=0.95, debug=FALSE, memory=NULL)
#' @examples
#' data = readseq(dataFile=getexample(dataType="protein"),dataType="protein")
#' data_process = pprocess(data=data,conservativeFilter=0.95)
#' @param data list, returns from the function "readseq()".
#' @param conservativeFilter numeric, a number in the range of 0~1,  0.95 by default. It's used to filter the highly conservative columns which the ratio of some residue is larger than the conservationFilter.
#' @param memory character, the type of matrix, NULL by default. It could be 'memory' or 'sparse'. If it's set to be 'memory', all data would be manipulated in the RAM by using normal matrix and package 'bigmemory'. If it's set to be 'sparse', the package "Matrix" would be used to manipulate massive matrices in memory and initialize huge sparse matrix, which could significantly reduce the RAM consumed. In default, it is set to be NULL, so that the system would determine automatically whether all data is manipulated in the RAM or not, according to the size of data inputted and the RAM available for R.
#' @param debug  logic, FALSE by default, indicates whether the debug message will be displayed or not.
#' @return list, contains the original data matrix, frequency matrix and other informations.
#' @references Du, X., Wang, Z., Wu, A., Song, L., Cao, Y., Hang, H., & Jiang, T. (2008). Networks of genomic co-occurrence capture characteristics of human influenza A (H3N2) evolution. Genome research, 18(1), 178-187. doi:10.1101/gr.6969007
#' @export
pprocess <- function(data=list(), conservativeFilter=0.95, debug=FALSE, memory=NULL){
  #parallel=FALSE
  cpus=NA
  t = Sys.time()
  cat("Pre-processing......")
  dataType = data$dataType
  conservativeFilter = cooccur.conservativeFilter(conservativeFilter,dataType,debug)

  threshold = conservativeFilter


  data <- cooccur.dataprepreprocess.preprocess(data=data, threshold=threshold, memory=memory, debug=debug )
  cat(sprintf("Pre-process......completed. (%s,%s)",nrow(data$matrix),ncol(data$matrix)))
  cooccur.printTimeCost('total co-occurrence.getco-occurrence time cost',t,debug)
  message("")
  return(data)
}

#' @title gencooccur
#'
#' @description Construct the co-occurrence network
#' @examples
#' data = readseq(dataFile=getexample(dataType="protein"),dataType="protein")
#' data_process = pprocess(data=data)
#' #cooccurNetwork  = gencooccur(data=data_process)
#' @param data list, returns from the function "pprocess()"
#' @param cooccurFilter numeric, a number in the range of 0~1. It determines whether two columns are perfect co-occurrence. In default, for the data type of protein, it is set to be 0.9, while for the other data types, it is set to be 1.
#' @param networkFile character, 'cooccurNetwork' be default. It is a file name with full path for storing the co-occurrence network for each row.
#' @param module logic, FALSE by default, to check whether the modules in each network of the networkFile would be calculated.
#' @param moduleFile character, 'cooccurNetworkModule' by default. It is a file name with full path for storing the modules for co-occurrence network.
#' @param property logic, FALSE by default, to check whether the properties for each network of the networkFile, including the network diameter, connectivity, ConnectionEffcient and so on, would be calculated.
#' @param propertyFile character, 'cooccurNetworkProperty' by default. It is a file name with full path storing the properties for each network of the networkFile.
#' @param siteCo  logic, FALSE by default, to check whether the residue co-occurence file would be calculated.
#' @param siteCoFile character, 'siteCooccurr' by default. It is a file name with full path for storing the RCOS between all pairs of columns, and the related p-values.
#' @param sampleTimes  numeric, an integer of permutations in the simulation when calculating the p-values. It should be greater than 100.
#' @param debug  logic, FALSE by default, indicates whether the debug message will be displayed or not.
#' @param parallel logic, FALSE by default. It only supports Unix/Mac (not Windows) system.
#' @return list, all the output file paths are attributed in it. \cr  The attribute "networkFile" stores the co-occurrence network for each row; \cr The attribute "moduleFile" is optional. When the module is set to be TRUE, it would be output. It stores the modules for co-occurrence network; \cr The attribute "propertyFile" is optional. When the property is set to be TRUE, it would be output. It stores the properties for co-occurrence network; \cr The attribute "siteCoFile" is optional. When the property is set to be TRUE, it would be output. It stores all the pairwise siteCos between columns.
#' @references Du, X., Wang, Z., Wu, A., Song, L., Cao, Y., Hang, H., & Jiang, T. (2008). Networks of genomic co-occurrence capture characteristics of human influenza A (H3N2) evolution. Genome research, 18(1), 178-187. doi:10.1101/gr.6969007
#' @export
gencooccur <- function(data=list(), cooccurFilter=NULL, networkFile='cooccurNetwork', module=FALSE,moduleFile='cooccurNetworkModule', property=FALSE, propertyFile='cooccurNetworkProperty', siteCo=FALSE, siteCoFile='siteCooccurr', sampleTimes=100, debug=FALSE, parallel=FALSE){
  t = Sys.time()
  #sequences <- cooccur.gennetework.calculateNetWork(sequences, alpha=alpha, parallel=parallel, cpus=cpus, filterfile=filterfile, rawfile=rawfile, modulefile=modulefile, propertyfile=propertyfile, cooccurfile=cooccurfile, pvaluefile=pvaluefile, ptimes=ptimes,debug=debug)
  #cooccur.printTimeCost('total co-occurrence.getco-occurrence time cost',t,debug)
  #return(sequences)
  #parallel=FALSE
  #cpus=NA
  dataType = data$dataType

  cooccurFilter = cooccur.cooccurFilter(cooccurFilter,dataType,debug)

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



  data <- cooccur.gennetework.calculateNetWork(sequences=data, alpha=alpha, parallel=parallel, filterfile=filterfile, rawfile=rawfile, modulefile=modulefile, propertyfile=propertyfile, cooccurfile=cooccurfile, pvaluefile=pvaluefile, ptimes=ptimes,debug=debug)
  cooccur.printTimeCost('total co-occurrence.createco-occurrence time cost',t,debug)
  message("")
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
#' @param dataFile character, a FASTA data file name with full path.
#' @param dataType character, 'protein' by default, the type of data will be processed. It could be 'DNA', 'RNA', 'protein', 'SNP' or 'other'.
#' @param conservativeFilter numeric, a number in the range of 0~1,  0.95 by default. It's used to filter the highly conservative columns which the ratio of some residue is larger than the conservationFilter.
#' @param cooccurFilter numeric, a number in the range of 0~1. It determines whether two columns are perfect co-occurrence. In default, for the data type of protein, it is set to be 0.9, while for the other data types, it is set to be 1.
#' @param networkFile character, 'cooccurNetwork' be default. It is a file name with full path for storing the co-occurrence network for each row.
#' @param module logic, FALSE by default, to check whether the modules in each network of the networkFile would be calculated.
#' @param moduleFile  character, 'cooccurNetworkModule' by default. It is a file name with full path for storing the modules for co-occurrence network.
#' @param property logic, FALSE by default, to check whether the properties for each network of the networkFile, including the network diameter, connectivity, ConnectionEffcient and so on, would be calculated.
#' @param propertyFile character, 'cooccurNetworkProperty' by default. It is a file name with full path storing the properties for each network of the networkFile.
#' @param siteCo logic, FALSE by default, to check whether the residue co-occurence file would be calculated.
#' @param siteCoFile character, 'siteCooccurr' by default. It is a file name with full path for storing the RCOS between all pairs of columns, and the related p-values.
#' @param sampleTimes  numeric, an integer of permutations in the simulation when calculating the p-values. It should be greater than 100.
# networkEvaluate    FALSE by default. If it is set to be TRUE, the statistical significance of all the edges in all the networks would be calculated by permutating the data as follows. Firstly, all columns in the data are randomly permutated; then, the co-occurrence network would be constructed. This process would be repeated N times (the value depends on the parameter sampleTimes). The p-value for each edge in original co-occurrence network was calculated as the number of times this edge appeared in the permutated networks. An additional file, named as "cooccurNetwork_pvalue" in default, would save the all the co-occurrence networks with p-value for all edges. Be attention! This process would be rather time-consuming.
#' @param parallel logic, FALSE by default. It only supports Unix/Mac (not Windows) system.
#' @param memory character, the type of matrix, NULL by default. It could be 'memory' or 'sparse'. If it's set to be 'memory', all data would be manipulated in the RAM by using normal matrix and package 'bigmemory'. If it's set to be 'sparse', the package "Matrix" would be used to manipulate massive matrices in memory and initialize huge sparse matrix, which could significantly reduce the RAM consumed. In default, it is set to be NULL, so that the system would determine automatically whether all data is manipulated in the RAM or not, according to the size of data inputted and the RAM available for R.
#' @param debug  logic, FALSE by default, indicates whether the debug message will be displayed or not.
#' @return list, all the output file paths are attributed in it. \cr  The attribute "networkFile" stores the co-occurrence network for each row; \cr The attribute "moduleFile" is optional. When the module is set to be TRUE, it would be output. It stores the modules for co-occurrence network; \cr The attribute "propertyFile" is optional. When the property is set to be TRUE, it would be output. It stores the properties for co-occurrence network; \cr The attribute "siteCoFile" is optional. When the property is set to be TRUE, it would be output. It stores all the pairwise siteCos between columns.
#' @references Du, X., Wang, Z., Wu, A., Song, L., Cao, Y., Hang, H., & Jiang, T. (2008). Networks of genomic co-occurrence capture characteristics of human influenza A (H3N2) evolution. Genome research, 18(1), 178-187. doi:10.1101/gr.6969007
#' @export
coocnet <-  function(dataFile="", dataType="protein", conservativeFilter=0.95, cooccurFilter=NULL, networkFile='cooccurNetwork', module=FALSE,moduleFile='cooccurNetworkModule', property=FALSE, propertyFile='cooccurNetworkProperty', siteCo=FALSE, siteCoFile='siteCooccurr', sampleTimes=100, debug=FALSE, parallel=FALSE, memory=NULL){
  #parallel=FALSE
  #cpus=NA
  t = Sys.time()
  cat(("reading file......"))
  sequences <- cooccur.readfasta.read(dataFile=dataFile, debug=debug)
  cat((sprintf("completed. (%s,%s)",nrow(sequences$matrix),ncol(sequences$matrix))))
  message("")
  cat("preprocess file......")

  conservativeFilter = cooccur.conservativeFilter(conservativeFilter,dataType,debug)
  #print(conservativeFilter)
  cooccurFilter = cooccur.cooccurFilter(cooccurFilter,dataType,debug)

  threshold = conservativeFilter
  #memory=NA
  sequences <- cooccur.dataprepreprocess.preprocess(data=sequences,threshold=threshold, memory=memory, debug=debug)
  cat(sprintf("preprocess file......completed. (%s,%s)",nrow(sequences$matrix),ncol(sequences$matrix)))
  message("")

  alpha = cooccurFilter
  ptimes = sampleTimes

  filterfile = NA
  cooccurfile=NA
  #cooccurfile="cooccurfile"

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

  networkEvaluate=FALSE
  if(networkEvaluate==TRUE){
    networkpfile = paste(networkFile, "_pvalue", sep="")
  }else{
    networkpfile = NA
  }

  if(siteCo==TRUE){
    pvaluefile = siteCoFile
  }else{
    pvaluefile=NA
  }

  sequences <- cooccur.gennetework.calculateNetWork(sequences=sequences, alpha=alpha, parallel=parallel, filterfile=filterfile, rawfile=rawfile, modulefile=modulefile, propertyfile=propertyfile, cooccurfile=cooccurfile, pvaluefile=pvaluefile, networkpfile=networkpfile, ptimes=ptimes,debug=debug)
  cooccur.printTimeCost('total co-occurrence.createco-occurrence time cost',t,debug)
  message("")
  return(sequences)


}




#' @title siteco
#'
#' @description Read and preprocess data, and calculate the pairwise site co-occurrence in one step
#' @examples
#' #pairwiseCooccur = siteco(dataFile=getexample(dataType="protein"), dataType="protein")
#' @param dataFile character, a FASTA data file name with full path.
#' @param dataType character, 'protein' by default, the type of data will be processed. It could be 'DNA', 'RNA', 'protein', 'SNP' or 'other'.
#' @param conservativeFilter numeric, a number in the range of 0~1,  0.95 by default. It's used to filter the highly conservative columns which the ratio of some residue is larger than the conservationFilter.
#' @param cooccurFilter numeric, a number in the range of 0~1. It determines whether two columns are perfect co-occurrence. In default, for the data type of protein, it is set to be 0.9, while for the other data types, it is set to be 1.
#' @param siteCoFile character, 'siteCooccurr' by default. It is a file name with full path for storing the RCOS between all pairs of columns, and the related p-values.
#' @param sampleTimes  numeric, an integer of permutations in the simulation when calculating the p-values. It should be greater than 100.
#' @param parallel logic, FALSE by default. It only supports Unix/Mac (not Windows) system.
#' @param memory character, the type of matrix, NULL by default. It could be 'memory' or 'sparse'. If it's set to be 'memory', all data would be manipulated in the RAM by using normal matrix and package 'bigmemory'. If it's set to be 'sparse', the package "Matrix" would be used to manipulate massive matrices in memory and initialize huge sparse matrix, which could significantly reduce the RAM consumed. In default, it is set to be NULL, so that the system would determine automatically whether all data is manipulated in the RAM or not, according to the size of data inputted and the RAM available for R.
#' @param debug  logic, FALSE by default, indicates whether the debug message will be displayed or not.
#' @return list, the output file path of 'sitecoFile' is attributed in it. The file stores all the pairwise siteCos between columns.
#' @references Du, X., Wang, Z., Wu, A., Song, L., Cao, Y., Hang, H., & Jiang, T. (2008). Networks of genomic co-occurrence capture characteristics of human influenza A (H3N2) evolution. Genome research, 18(1), 178-187. doi:10.1101/gr.6969007
#' @export
siteco <-  function(dataFile="", dataType="protein", conservativeFilter=0.95, cooccurFilter=NULL, siteCoFile='siteCooccurr', sampleTimes=100, debug=FALSE, parallel=FALSE, memory=NULL){
  #parallel=FALSE
  cpus=NA
  t = Sys.time()
  cat(("reading file......"))
  sequences <- cooccur.readfasta.read(dataFile=dataFile, debug=debug)
  cat((sprintf("completed. (%s,%s)",nrow(sequences$matrix),ncol(sequences$matrix))))
  message("")
  cat("preprocess file......")

  conservativeFilter = cooccur.conservativeFilter(conservativeFilter,dataType,debug)

  cooccurFilter = cooccur.cooccurFilter(cooccurFilter,dataType,debug)

  threshold = conservativeFilter

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

  sequences <- cooccur.gennetework.calculateNetWork(sequences=sequences, alpha=alpha, parallel=parallel, filterfile=filterfile, rawfile=rawfile, modulefile=modulefile, propertyfile=propertyfile, cooccurfile=cooccurfile, pvaluefile=pvaluefile, ptimes=ptimes,debug=debug)
  cooccur.printTimeCost('total co-occurrence.createco-occurrence time cost',t,debug)
  message("")
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
#' @param networkFile character, a file name with full path for storing the co-occurrence network for each row and generated from cooc() or gencooccur().
#' @param networkNames character, a vector of network names
#' @return   a list of igraph::graph.data.frame object ordered by the input network names.
#' @export
toigraph <-  function(networkFile="", networkNames=c()){

  if(!requireNamespace("igraph", quietly = TRUE)){
    stop("Package 'igraph' is required.")
  }

  dataframe = read.csv(networkFile,sep=" ",header=FALSE,stringsAsFactors = FALSE)
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
  message("")
  return(graphList)
}





#'getexample
#' @title getexample
#' @description Get the example data
#' @note  Both the dataset "DNA" and "protein" are sampled from the Hemagglutinin sequences of human influenza H3N2 viruses, while the dataset "SNP" are simulated.
#' @usage getexample(dataType)
#' @examples
#' dataFile = getexample(dataType='protein')
#' @param dataType character, 'protein' by default. It could be 'DNA', 'RNA', 'protein', 'SNP' or 'other'.
#' @return character, the file path of the example data
#' @export
#'
getexample <- function(dataType='protein'){
  path = ""
  if(dataType=="DNA"){
    path = system.file("extdata", "DNA", package = "cooccurNet")
  }else if(dataType=="protein"){
    path = system.file("extdata", "protein", package = "cooccurNet")
  }else if(dataType=="SNP"){
    path = system.file("extdata", "SNP", package = "cooccurNet")
  }else if(dataType=="RNA"){
    path = system.file("extdata", "RNA", package = "cooccurNet")
  }else{
    stop("the param dataType should in the range of [DNA,RNA,protein,SNP]")
  }
  message("")
  return(path)
}


#'changeLog
#' @title changeLog
#' @description Get the most recent n lines
#' @note  Once you have installed cooccurNet, the change-log can also be viewed from the R prompt.
#' @usage changeLog(n=20)
#' @examples
#' logs = changeLog(n=20)
#' @param n numeric, 20 by default, the number of lines will be shown up.
#' @return list, the most recent n lines
#' @export
#'
changeLog <- function(n=20){
  path = ""
  path = system.file("extdata", "changeLog", package = "cooccurNet")
  fileName<-paste(path, "changeLog" ,sep="/")
  #print(path)
  logs = readLines(fileName)
  if(n>length(logs)){
    n = length(logs)
  }
  message("")
  return(logs[1:n])
}



#'getexample_forRCOS
#' @title getexample_forRCOS
#' @description Get the test data for testing RCOS method
#' @note  File paths for all the files available for testing RCOS. Currently, there is only one file "HA1protein_humanH3N2", the sequences within which are derived from the database of Influenza Virus Resource.
#' @usage getexample_forRCOS()
#' @examples
#' dataFile = getexample_forRCOS()
#' @return character, the files available for testing RCOS.
#' @export
#'
getexample_forRCOS <- function(){
  path = ""
  path = system.file("extdata", "RCOS", package = "cooccurNet")
  fileName<-paste(path, dir(path) ,sep="/")
  #print(path)
  message("")
  return(fileName)
}
