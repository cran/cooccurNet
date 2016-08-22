\name{cooccurNet}
\alias{cooccurNet}
\title{cooccurNet: Co-Occurrence network computation for R}

\description{
Read and preprocess fasta format data, and construct the co-occurrence network for downstream analyses. This R package is to construct the co-occurrence network with the algorithm developed by Du (Du 2008). It could be used to transform the data with high-dimension, such as DNA or protein sequence, into co-occurrence networks. Co-occurrence network could not only capture the co-variation pattern between variables, such as the positions in DNA or protein sequences, but also reflect the relationship between samples. Although it is originally used in DNA and protein sequences, it could be also used to other kinds of data, such as SNP.
}

\details{
Index of functions/methods (grouped in a friendly way):

1. readseq(dataFile="",  dataType="DNA", debug=FALSE)

2. pprocess(data=list(),  conservativeFilter=0.95, debug=FALSE)

3. gencooccur(data=list(), cooccurFilter, networkFile='cooccurNetwork', \cr module=FALSE,moduleFile='cooccurNetworkModule', property=FALSE, propertyFile='cooccurNetworkProperty', siteCo=FALSE, siteCoFile='siteCooccurr', sampleTimes=100, debug=FALSE)

4. coocnet(dataFile, dataType="DNA", conservativeFilter=NULL, cooccurFilter=NULL, networkFile='cooccurNetwork', module=FALSE,moduleFile='cooccurNetworkModule', property=FALSE, propertyFile='cooccurNetworkProperty', siteCo=FALSE, siteCoFile='siteCooccurr', sampleTimes=100,  debug=FALSE)

5. siteco(dataFile, dataType="DNA", conservativeFilter=NULL, cooccurFilter=NULL, siteCoFile='siteCooccurr', sampleTimes=100, debug=FALSE)

6. toigraph(networkFile="", networkNames=c())

7. getexample(dataType)

}

\note{

Arguments:

1. data \cr
type: list \cr
description: it contains the DNA, protein, SNP data or other kinds of data read from the input data file.


2. dataFile \cr
type: character \cr
description: the file with full path of the DNA, protein, SNP data or other kinds of data.


3. conservativeFilter \cr
type: numeric \cr
description: a value in the range of 0~1. The column with conservative score greater than it would be filtered in the later analyses.

4. cooccurFilter \cr
type: numeric \cr
description: a value in the range of 0~1. It determines whether two columns are perfect co-occurrence.

5. network \cr
type: logical \cr
description: The default value is TRUE. It determines whether the network file is to be generated.

6. networkFile \cr
type: character \cr
description: the file with full path for storing the co-occurrence network for each row;

7. module \cr
type: logical \cr
description: The default value is FALSE. If it is set to be TRUE, the modules in each network of the networkFile would be calculated.

8. moduleFile \cr
type: character \cr
description: the file with full path for storing the modules for co-occurrence network;

9. property \cr
type: logical \cr
description: The default value is FLASE. If it is set to be TRUE, the properties for each network of the networkFile, including the network diameter, connectivity, ConnectionEffcient and so on, would be calculated.

10. propertyFile \cr
type: character \cr
description: the file with full path for storing the modules for co-occurrence network;

11. siteCo \cr
type: logical \cr
description: The default value is FALSE. If it is set to be TRUE, the extent of co-occurrence between all pairs of columns would be calculated. It is defined as the ratio of rows with perfect co-occurrence.

12. siteCoFile \cr
type: character \cr
description: file for storing the extent of co-occurrence between all pairs of columns, and the related p-values. The later are calculated by simulations as follows: firstly, all columns in the data are randomly permutated; then, the pairwise siteCos are calculated. This process would be repeated N times (the value depends on the parameter sampleTimes). For each pair of columns, the rank of the original siteCo in the N siteCos derived from simulations are considered as the p-value for the original siteCo.

13. sampleTimes \cr
type: integer \cr
description: the number of permutations in the simulation when calculating the p-values.

14.debug \cr
type: logical \cr
description: The default value is FALSE. It determines whether the debug messages is shown.


Output files:

1. networkFile \cr
filePath = data$networkFile \cr

2. moduleFile \cr
filePath = data$moduleFile \cr

3. propertyFile \cr
filePath = data$propertyFile \cr

4. siteCoFile \cr
filePath = data$siteCoFile \cr

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

  #example of readseq()
  #read sequences from the sample fasta file
  #data = readseq(dataFile=dataFile, dataType="protein")

  #example of pprocess()
  #preprocess the sequence dataFile
  #data_process = pprocess(data=data)

  #example of gencooccur()
  #generate co-occurrence network and save it into the 'networkFile'
  #cooccurNetwork  = gencooccur(data=data_process)
  #check the 'networkFile' path
  #print(cooccurNetwork$networkFile)

  #example of coocnet()
  #also, you can generate the co-occurrence network by one command
  cooccurNetwork  = coocnet(dataFile=dataFile, dataType="protein")
  #check the 'networkFile' path
  #print(cooccurNetwork$networkFile)

  #example of siteco()
  #you can generate the co-occurrence network siteCoFile by one command
  #this command will take long time to calculate the p-value.
  #pairwiseCooccur = siteco(dataFile=dataFile, dataType="protein")
  #check the 'siteCoFile' path
  #print(pairwiseCooccur$siteCoFile)

  #example of toigraph()
  #you can transform a network file to the igraph.data.frame
  #cooccurNetwork  = coocnet(dataFile=dataFile,dataType="protein")
  #get igraph data frame by specifying the network name
  #network_igraph = toigraph(networkFile=cooccurNetwork$networkFile, networkName=c("EPI823725"))
  #Plot the network (The package "igraph" must be installed and loaded firstly)
  #read the names of network
  #networkName = cooccurNetwork$xnames
  #Transform all cooccurrence network into the igraph format
  #Network_igraph = toigraph(networkFile=cooccurNetwork$networkFile, networkNames=networkName)



}
