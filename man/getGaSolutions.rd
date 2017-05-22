\name{getGaSolutions}
\alias{getGaSolutions}
\alias{Kmatfunc}
\alias{calculatecrossvalue} 
\alias{getstats}
\alias{mapfunct}
\title{getGaSolutions}
\description{
 Selection in breeding programs can be done by using phenotypes (phenotypic selection), pedigree relationship (breeding value selection) or molecular markers (marker assisted selection or genomic selection). All these methods are based on truncation selection, focusing on the best performance of parents before mating. In this article we proposed an approach to breeding, named genomic mating, which focuses on mating instead of truncation selection. Genomic mating uses information in a similar fashion to genomic selection but includes information on complementation of parents to be mated. Following the efficiency frontier surface, genomic mating uses concepts of estimated breeding values, risk (usefulness) and coefficient of ancestry to optimize mating between parents. This package uses a genetic algorithm to find solutions to this optimization problem and the results from our simulations comparing genomic selection, phenotypic selection and the mating approach indicate that current approach for breeding complex traits is more favorable than phenotypic and genomic selection. Genomic mating is similar to genomic selection in terms of estimating marker effects, but in genomic mating the genetic information and the estimated marker effects are used to decide which genotypes should be crossed to obtain the next breeding population.
}
\usage{
getGaSolutions(Markers, K, markereffects,nmates=NULL,minparents=10, 
impinbreedstepsize=.2, impvar, impforinbreed,npopGA, nitGA, 
plotiters,nelite, mutprob, mc.cores, miniters=100,
minitbefstop=80, tolparconv=1e-6)
}
\arguments{
  \item{Markers}{The matrix of markers rows corresponding to individuals and columns for markers, the markers scores are coded as -1,0,1.}
  \item{K}{symmetric genomic relationship matrix, the order of the row and columns of this matrix should follow the order of genotypes in the rows of \code{Markers}.}
  \item{markereffects}{effects of markers for a trait}
  \item{nmates}{number of mates to select, default value is NULL (number of mates is equal to number of mates)}
  \item{minparents}{minimum number of parents in the solution (importance parameter for inbreeding is increased till minimum number of parents are included in the mating solution), minimum is 1.}
  \item{impinbreedstepsize}{stepsize for importance parameter for inbreeding to be increased till minimum number of parents are included in the mating solution,}
  \item{impvar}{importance parameter of the cross variance term}
   \item{impforinbreed}{importance parameter for inbreeding}
  \item{npopGA}{genetic algorithm parameter: number of solutions generated at each cycle of the GA}
  \item{nitGA}{genetic algorithm parameter: number of GA cycles before algorithm stops }
  \item{plotiters}{genetic algorithm parameter: if TRUE the value of the objective function over iterations will be plotted}
  \item{nelite}{genetic algorithm parameter:number of elite solutions selected at each cycle of the GA}
  \item{mutprob}{genetic algorithm parameter: mutation probability}
  \item{mc.cores}{genetic algorithm parameter: number of cores to use}
  \item{miniters}{genetic algorithm parameter: minimum number of GA cycles before algorithm stops }
  \item{minitbefstop}{genetic algorithm parameter: minimum number of GA cycles before algorithm continues when the tolerance is reached (no change in the criterion value)}
  \item{tolparconv}{genetic algorithm parameter: the maximum change in criterion value accepted for conbergence.}

}
\value{
  Returns a list with two elements: the first element in this list is the list of mates in the best solution (in terms of row numbers corresponding to the genotypes in \code{Markers}), the second element is the criterion value for this solution.    
}
\details{
 The efficient mating problem can be stated as an optimization problem as follows: 
 minimize w.r.t. \eqn{P_{32}} \eqn{ r( \lambda_1, \lambda_2, P_{32}) =  -Risk(\lambda_1, P_{32}) + \lambda_2 * Inbreeding(P_{32})} where \eqn{\lambda_2\geq 0} is the parameter whose magnitude controls the amount of co-ancestry in the progeny, and the minimization is over the space of the mating matrices \eqn{P_{32}} construction of which is described  in detail below. \eqn{\lambda_1} controls allele heterozygosity weighted by the marker effects and  \eqn{\lambda_2} controls allele diversity. When \eqn{\lambda_1=0} the risk measure is the same as total expected gain. 



 }
\references{
 Akdemir, Deniz, and Julio I. Sanchez. ''Efficient Breeding by Genomic Mating.'' Frontiers in Genetics 7 (2016).
 
 }
\examples{
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
gasols<-getGaSolutions(Markers, K, markereffects,minparents=1, 
impinbreedstepsize=.2, impvar=.001, 
impforinbreed=10,npopGA=40, nitGA=10, plotiters=TRUE,
 mc.cores=1,nelite=5, mutprob=0.8)
gasols
#increase npopGA, nitGA for more reasonable solutions 
#(npopGA=40, nitGA=10) too low for any convergence
}
\author{Deniz Akdemir}

