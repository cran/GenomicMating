\name{Amat.pieces}
\alias{Amat.pieces}
\title{Amat.pieces}
\description{
This calculates the genomic relationship matrix using the formula in Van Ragen (2008)}
\usage{
Amat.pieces(M, pieces=10, mc.cores=1)
}
\arguments{
  \item{M}{The matrix of markers rows rorresponding to individuals and columns for markers, the markers scores are coded as 0,1,2.}
  \item{pieces}{number of chunks to split the markers}
  \item{mc.cores}{number of cores to use}
  }
\value{
 a genomic relationship matrix
}

\references{
VanRaden, Paul M. ''Efficient methods to compute genomic predictions.'' Journal of dairy science 91.11 (2008): 4414-4423.
 
 }
\examples{
library(GenomicMating)
N=50
nmarkers=500
Markers<-c()
for (i in 1:N){
  Markers<-rbind(Markers,sample(-1:1,nmarkers, replace=TRUE))
}

markereffects<-rep(0,nmarkers)
markereffects[sample(1:nmarkers,nmarkers/2)]<-rnorm(nmarkers/2)
Markers[1:5,1:5]

K=Amat.pieces(Markers+1, pieces=5) 
K[1:5,1:5]

}
\author{Deniz Akdemir}
