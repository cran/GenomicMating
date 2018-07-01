\name{getGaSolutionsFrontier}
\alias{getGaSolutionsFrontier}
\title{getGaSolutionsFrontier}
\description{Genomic mating is similar to genomic selection in terms of estimating marker effects, but in genomic mating the genetic information and the estimated marker effects are used to decide which genotypes should be crossed to obtain the next breeding population. This program uses genetic algorithm to obtain the frontier solutions that minimize inbreeding and maximize gain and usefulness. The solutions in the optimized frontier are nondominated.}
\usage{getGaSolutionsFrontier(Markers, Markers2=NULL,K, 
markereffects,markermap=NULL,nmates=NULL,npopGA=100, 
nitGA=100, mutprob=1, mc.cores=1, noself=F,method=1L, 
type=0L, generation=1L,plotiters=F)
}
\arguments{
\item{Markers}{The matrix of markers rows corresponding to individuals and columns for markers, the markers scores are coded as -1,0,1. (For Method=3 the markers are coded as probabilities between 0 and 1.)}
 \item{Markers2}{The matrix of markers rows corresponding to individuals and columns for markers, the markers scores are coded as -1,0,1. (For Method=3 the markers are coded as probabilities between 0 and 1.)}
\item{K}{symmetric genomic relationship matrix, the order of the row and columns of this matrix should follow the order of genotypes in the rows of \code{rbind(Markers, Markers2)}.}
\item{markermap}{A map for markers. two columns, firts column is named chr second named pos for the chromosome and position of the markers specified above}
\item{markereffects}{effects of markers for a trait}
\item{nmates}{number of mates to select, default value is NULL (number of mates is equal to number of mates)}
\item{npopGA}{genetic algorithm parameter: number of solutions generated at each cycle of the GA}
\item{nitGA}{genetic algorithm parameter: number of GA cycles before algorithm stops }
\item{mutprob}{genetic algorithm parameter: mutation probability}
\item{mc.cores}{genetic algorithm parameter: number of cores to use}
\item{plotiters}{genetic algorithm parameter: if TRUE the value of the objective function over iterations will be plotted}
\item{noself}{ Is selfing allowed? (TRUE or FALSE)}
\item{method}{Which method to use? (1,2,3) See Detalis.}
\item{type}{Only for method=2. Type of crosses (1 (DH), 2 (RISELF)).}
\item{generation}{Only for method=2. Generation at which the cross variances are calculated.}
}
\value{Returns a list with two elements: the first element in this list is a list of solutions found on the frontier, the second element is the matrix of criterion values (Gain, Usefulness, and Inbreeding) corresponding to these solutions.
}
\details{This program uses genetic algorithm to produce a number of solutions on the frontier curve simultaneously for the multi-objective optimization problem which is defined by minimization of \eqn{-Gain(P_{32}),} \eqn{-Usefulness(P_{32})} and \eqn{Inbreeding(P_{32})} with respect to \eqn{P_{32}.}
}
\references{
Akdemir,Sanchez. "Efficient Breeding by Genomic Mating." Frontiers in Genetics (2016).

Lehermeier at al. "Genetic gain increases by applying the usefulness criterion with improved variance prediction in selection of crosses" Genetics (2017).

Broman et al. "R/qtl: QTL mapping in experimental crosses." Bioinformatics (2003).

VanRaden, Paul M. ''Efficient methods to compute genomic predictions.'' Journal of dairy science (2008).
}
\examples{
	\dontrun{
library(GenomicMating)
#############
####for method 3 polyploid. Markers need to be coded between 0 and 1.
N=20
nmarkers=100
Markers<-c()
for (i in 1:N){
  Markers<-rbind(Markers,runif(nmarkers))
}


Markers2<-c()
for (i in 1:N){
  Markers2<-rbind(Markers2,runif(nmarkers))
}

markereffects<-rep(0,nmarkers)
markereffects[sample(1:nmarkers,nmarkers/2)]<-rnorm(nmarkers/2)
Markers[1:5,1:5]
#library(parallel)
K=Amat.pieces(rbind(Markers)*2-1, pieces=5) 

K2=Amat.pieces(rbind(Markers,Markers2)*2-1, pieces=5) 
K[1:5,1:5]
rownames(Markers)<-paste("l", 1:nrow(Markers),sep="_")
rownames(Markers2)<-paste("l", (nrow(Markers)+1):(nrow(Markers)+
nrow(Markers2)),sep="_")
rownames(K2)<-colnames(K2)<-c(rownames(Markers),rownames(Markers2))
rownames(K)<-colnames(K)<-c(rownames(Markers))

which.max(Markers\%*\%markereffects)
markermap=as.matrix(data.frame(chr=rep(1,nmarkers),
pos=seq(0,1,length=nmarkers)))

colnames(Markers)<-1:nmarkers


gasols4<-getGaSolutionsFrontier(Markers=Markers,Markers2=Markers2, K=K2,
markereffects,markermap=markermap,nmates=10,npopGA=100, nitGA=100,
mc.cores=1, mutprob=0.999, noself= TRUE, method=3,
type=2L, generation=1L, plotiters= TRUE)



###plot results

pairs(gasols4[[1]])

####Use plotGM.

plotGM(GMsols=gasols4, type="3D", traitnum=1)
plotGM(GMsols=gasols4, type="SOM", traitnum=1)


}
}
\author{Deniz Akdemir, Julio Isidro Sanch\'ez, Hanna Haikka, Itaraju Baracuhy Brum}

