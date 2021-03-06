% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cooccurrence.R
\name{getexample}
\alias{getexample}
\title{getexample}
\usage{
getexample(dataType)
}
\arguments{
\item{dataType}{character, 'protein' by default. It could be 'DNA', 'RNA', 'protein', 'SNP' or 'other'.}
}
\value{
character, the file path of the example data
}
\description{
Get the example data
}
\details{
getexample
}
\note{
Both the dataset "DNA" and "protein" are sampled from the Hemagglutinin sequences of human influenza H3N2 viruses, while the dataset "SNP" are simulated.
}
\examples{
dataFile = getexample(dataType='protein')
}

