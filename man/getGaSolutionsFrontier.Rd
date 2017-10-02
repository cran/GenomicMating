\name{getGaSolutionsFrontier}
\alias{getGaSolutionsFrontier}
\title{getGaSolutionsFrontier}
\description{
 Selection in breeding programs can be done by using phenotypes (phenotypic selection), pedigree relationship (breeding value selection) or molecular markers (marker assisted selection or genomic selection). All these methods are based on truncation selection, focusing on the best performance of parents before mating. In this article we proposed an approach to breeding, named genomic mating, which focuses on mating instead of truncation selection. Genomic mating uses information in a similar fashion to genomic selection but includes information on complementation of parents to be mated. Following the efficiency frontier surface, genomic mating uses concepts of estimated breeding values, risk (usefulness) and coefficient of ancestry to optimize mating between parents. This package uses a genetic algorithm to find solutions to this optimization problem and the results from our simulations comparing genomic selection, phenotypic selection and the mating approach indicate that current approach for breeding complex traits is more favorable than phenotypic and genomic selection. Genomic mating is similar to genomic selection in terms of estimating marker effects, but in genomic mating the genetic information and the estimated marker effects are used to decide which genotypes should be crossed to obtain the next breeding population.
}
\usage{
getGaSolutionsFrontier(Markers, K, markereffects,nmates=NULL,
npopGA, nitGA, mutprob, mc.cores)
}
\arguments{
  \item{Markers}{The matrix of markers rows corresponding to individuals and columns for markers, the markers scores are coded as -1,0,1.}
  \item{K}{symmetric genomic relationship matrix, the order of the row and columns of this matrix should follow the order of genotypes in the rows of \code{Markers}.}
  \item{markereffects}{effects of markers for a trait}
  \item{nmates}{number of mates to select, default value is NULL (number of mates is equal to number of mates)}
  \item{npopGA}{genetic algorithm parameter: number of solutions generated at each cycle of the GA}
  \item{nitGA}{genetic algorithm parameter: number of GA cycles before algorithm stops }
  \item{mutprob}{genetic algorithm parameter: mutation probability}
  \item{mc.cores}{genetic algorithm parameter: number of cores to use}
  }
\value{
  Returns a list with two elements: the first element in this list is a list of solutions found on the frontier, the second element is the matrix of criterion values (Gain, Usefulness, and Inbreeding) corresponding to these solutions.    
}
\details{
 This program uses genetic algorithm to  produce a number of solutions on the frontier curve simultaneously for the multi-objective optimization problem which is defined by minimization of \eqn{-Gain(P_{32}),} \eqn{-Usefulness(P_{32})} and \eqn{Inbreeding(P_{32})} with respect to \eqn{P_{32}.}



 }
\references{
 Akdemir, Deniz, and Julio I. Sanchez. ''Efficient Breeding by Genomic Mating.'' Frontiers in Genetics 7 (2016).
 
 }
\examples{
	\dontrun{
library(GenomicMating)
N=50
nmarkers=1000
Markers<-c()
for (i in 1:N){
  Markers<-rbind(Markers,sample(-1:1,nmarkers, replace=TRUE))
}

markereffects<-rep(0,nmarkers)
markereffects[sample(1:nmarkers,nmarkers/2)]<-rnorm(nmarkers/2)
Markers[1:5,1:5]

K=Amat.pieces(Markers+1, pieces=5) 
K[1:5,1:5]

rownames(Markers)<-rownames(K)<-colnames(K)<-paste("l", 1:nrow(K),sep="_")

gasols2<-getGaSolutionsFrontier(Markers, K, markereffects,
				nmates=10,npopGA=10, nitGA=40,
                 mc.cores=1, mutprob=0.999)

str(gasols2)
library(scatterplot3d)
scatterplot3d(gasols2[[1]])

#increase npopGA, nitGA for more reasonable solutions 
#(npopGA=40, nitGA=10) too low for any convergence
}
}
\author{Deniz Akdemir}

