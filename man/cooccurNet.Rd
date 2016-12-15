\name{cooccurNet}
\alias{cooccurNet}
\title{cooccurNet: Co-Occurrence network computation for R}

\description{
Read and preprocess fasta format data, and construct the co-occurrence network for downstream analyses. This R package is to construct the co-occurrence network with the algorithm developed by Du (2008) <DOI:10.1101/gr.6969007>. It could be used to transform the data with high-dimension, such as DNA or protein sequence, into co-occurrence networks. Co-occurrence network could not only capture the co-variation pattern between variables, such as the positions in DNA or protein sequences, but also reflect the relationship between samples. Although it is originally used in DNA and protein sequences, it could be also used to other kinds of data, such as RNA, SNP, etc.
}

\details{
Index of functions/methods (grouped in a friendly way):

1. readseq(dataFile="",  dataType="protein", debug=FALSE)

2. pprocess(data=list(), conservativeFilter=0.95, memory=NULL, debug=FALSE)

3. gencooccur(data=list(), cooccurFilter=NULL, networkFile='cooccurNetwork', module=FALSE, moduleFile='cooccurNetworkModule', property=FALSE, propertyFile='cooccurNetworkProperty', siteCo=FALSE, siteCoFile='siteCooccurr', sampleTimes=100, debug=FALSE, parallel=FALSE)

4. coocnet(dataFile="", dataType="protein", conservativeFilter=0.95, memory=NULL, cooccurFilter=NULL, networkFile='cooccurNetwork', module=FALSE,moduleFile='cooccurNetworkModule', property=FALSE, propertyFile='cooccurNetworkProperty', siteCo=FALSE, siteCoFile='siteCooccurr', sampleTimes=100, debug=FALSE, parallel=FALSE)

5. siteco(dataFile="", dataType="protein", conservativeFilter=0.95, memory=NULL, cooccurFilter=NULL, siteCoFile='siteCooccurr', sampleTimes=100, debug=FALSE, parallel=FALSE)

6. toigraph(networkFile="", networkNames=c())

7. getexample(dataType)

8. getexample_forRCOS(dataType)

9. changeLog(n=20)

}

\note{

Arguments:

1. dataFile \cr
type: character \cr
description: a FASTA data file name with full path.

2. dataType \cr
type: character \cr
description: 'protein' by default, the type of data will be processed. It could be 'DNA', 'RNA', 'protein', 'SNP' or 'other'.

3. conservativeFilter \cr
type: numeric \cr
description:  a number in the range of 0~1,  0.95 by default. It's used to filter the highly conservative columns which the ratio of some residue is larger than the conservationFilter.

4. cooccurFilter \cr
type: numeric \cr
description:  a number in the range of 0~1. It determines whether two columns are perfect co-occurrence. In default, for the data type of protein, it is set to be 0.9, while for the other data types, it is set to be 1.

5. networkFile \cr
type: character \cr
description: 'cooccurNetwork' be default. It is a file name with full path for storing the co-occurrence network for each row.

6. module \cr
type: logical \cr
description: FALSE by default, to check whether the modules in each network of the networkFile would be calculated.

7. moduleFile \cr
type: character \cr
description:  'cooccurNetworkModule' by default. It is a file name with full path for storing the modules for co-occurrence network.

8. property \cr
type: logical \cr
description:  FALSE by default, to check whether the properties for each network of the networkFile, including the network diameter, connectivity, ConnectionEffcient and so on, would be calculated.

9. propertyFile \cr
type: character \cr
description:  'cooccurNetworkProperty' by default. It is a file name with full path storing the properties for each network of the networkFile.

10. siteCo \cr
type: logical \cr
description: FALSE by default, to check whether the residue co-occurence file would be calculated.

11. siteCoFile \cr
type: character \cr
description:  'siteCooccurr' by default. It is a file name with full path for storing the RCOS between all pairs of columns, and the related p-values.

12. sampleTimes \cr
type: numeric \cr
description: a integer of permutations in the simulation when calculating the p-values. It should be greater than 100.

13.debug \cr
type: logical \cr
description: FALSE by default, indicates whether the debug message will be displayed or not.

14.memory \cr
type: character \cr
description:  the type of matrix, NULL by default. It could be 'memory' or 'sparse'. If it's set to be 'memory', all data would be manipulated in the RAM by using normal matrix and package 'bigmemory'. If it's set to be 'sparse', the package "Matrix" would be used to manipulate massive matrices in memory and initialize huge sparse matrix, which could significantly reduce the RAM consumed. In default, it is set to be NULL, so that the system would determine automatically whether all data is manipulated in the RAM or not, according to the size of data inputted and the RAM available for R.

15.parallel \cr
type: logical \cr
description: FALSE by default. It only supports Unix/Mac (not Windows) system.

Output files:

1. networkFile \cr
filePath = data$networkFile \cr
description: A co-occurrence network (networkFile) was given for each sequence, which could be easily transformed into the format of igraph by the function toigraph().

2. moduleFile \cr
filePath = data$moduleFile \cr
description: The inherent modules(moduleFile) in each co-occurrence network were given for each sequence.

3. propertyFile \cr
filePath = data$propertyFile \cr
description: Some basic network attributes(propertyFile) such as connectivity and clustering coefficient for each network were given.

4. siteCoFile \cr
filePath = data$siteCoFile \cr
description: The extent of co-occurrence between residues(siteCoFile), defined as residue co-occurrence score (RCOS), were given for all pairs of residues.

For more details, please see the 'Vignette: Extending cooccurNet' by using the following command.
vignette("Extending-cooccurNet",package="cooccurNet")

}
\author{

Yuanqiang Zou, Yousong Peng, and Taijiao Jiang

Maintainers: Yuanqiang Zou <jerrytsou2001@gmail.com>

}

\references{
   Du, X., Wang, Z., Wu, A., Song, L., Cao, Y., Hang, H., & Jiang, T. (2008). Networks of genomic co-occurrence capture characteristics of human influenza A (H3N2) evolution. Genome research, 18(1), 178-187. doi:10.1101/gr.6969007
}

\examples{
  #example of get example data
  dataFile=getexample(dataType="protein")

  #example of get file paths for all the files available for testing RCOS.
  #dataFiles = getexample_forRCOS()

  #example of readseq()
  #read sequences from the sample fasta file
  #data = readseq(dataFile=dataFile, dataType="protein")

  #example of pprocess()
  #preprocess the sequence dataFile
  #data_process = pprocess(data=data,conservativeFilter=0.95)

  #example of gencooccur()
  #generate co-occurrence network and save it into the 'networkFile'
  #cooccurNetwork  = gencooccur(data=data_process, cooccurFilter=0.9, networkFile='cooccurNetwork')
  #check the 'networkFile' path
  #print(cooccurNetwork$networkFile)

  #example of coocnet()
  #also, you can generate the co-occurrence network by one command
  #cooccurNetwork  = coocnet(dataFile=getexample(dataType="protein"), dataType="protein")
  #check the 'networkFile' path
  #print(cooccurNetwork$networkFile)

  #example of siteco()
  #you can generate the co-occurrence network siteCoFile by one command
  #this command will take long time to calculate the p-value.
  #pairwiseCooccur = siteco(dataFile=getexample(dataType="protein"), dataType="protein")
  #check the 'siteCoFile' path
  #print(pairwiseCooccur$siteCoFile)

  #example of toigraph()
  #you can transform a network file to the igraph.data.frame
  #cooccurNetwork  = coocnet(dataFile=getexample(dataType="protein"),dataType="protein")
  #get igraph data frame by specifying the network name
  #network_igraph = toigraph(networkFile=cooccurNetwork$networkFile, networkName=c("EPI823725"))
  #Plot the network (The package "igraph" must be installed and loaded firstly)
  #read the names of network
  #networkName = cooccurNetwork$xnames
  #Transform all cooccurrence network into the igraph format
  #Network_igraph = toigraph(networkFile=cooccurNetwork$networkFile, networkNames=networkName)

  #example of changelog()
  #logs = changeLog(n=20)


}
