\name{getGaSolutionsFrontierMultiTraitSimcross}
\alias{getGaSolutionsFrontierMultiTraitSimcross}
\title{getGaSolutionsFrontierMultiTraitSimcross}
\description{A generalization of the single trait mating approach in getGaSolutionsFrontier. The usefulness statistic is calculated using the simcross function in the qtl package.}
\usage{getGaSolutionsFrontierMultiTraitSimcross(Markers, Markers2=NULL,K,
map,markereffectslist,nmates=NULL,  nSim = 5,npopGA, nitGA,
mutprob, mc.cores, noself=F,simtype="riself", plotiters=F)}
\arguments{\item{Markers}{The matrix of markers rows corresponding to individuals and columns for markers, the markers scores are coded as -1,0,1. (For Method=3 the markers are coded as probabilities between 0 and 1.)}
 \item{Markers2}{The matrix of markers rows corresponding to individuals and columns for markers, the markers scores are coded as -1,0,1. (For Method=3 the markers are coded as probabilities between 0 and 1.)}
  \item{K}{symmetric genomic relationship matrix, the order of the row and columns of this matrix should follow the order of genotypes in the rows of \code{rbind(Markers, Markers2)}.}
 \item{map}{a map for markers. two columns, first column is named chr second named pos for the chromosome and position of the markers specified above.}
 \item{markereffectslist}{effects of markers for several traits given as a list.}
  \item{nmates}{number of mates to select, default value is NULL (number of mates is equal to number of genotypes).}
   \item{nSim}{number of progeny simulated for each pair.}
  \item{npopGA}{genetic algorithm parameter: number of solutions generated at each cycle of the GA.}
  \item{nitGA}{genetic algorithm parameter: number of GA cycles before algorithm stops.}
  \item{mutprob}{genetic algorithm parameter: mutation probability.}
  \item{mc.cores}{genetic algorithm parameter: number of cores to use.}
  \item{noself}{Default is FALSE. Specifies if selfing is allowed.}
  \item{simtype}{Default is ''riself''. Argument passed to simcross, not all types work with the GenomicMating package.}
  \item{plotiters}{Logical. Default is FALSE. Iterations are plotted if TRUE.}
  }
\value{Returns a list with two elements: the first element in this list is a list of solutions found on the frontier, the second element is the matrix of criterion values (Gain, Usefulness, and Inbreeding) corresponding to these solutions.}
\details{This program uses genetic algorithm to  produce a number of solutions on the frontier curve simultaneously for the multi-objective optimization problem which is defined by minimization of \eqn{-Gain(P_{32}),} \eqn{-Usefulness(P_{32})} and \eqn{Inbreeding(P_{32})} with respect to \eqn{P_{32}.}}
\references{Akdemir,Sanchez. "Efficient Breeding by Genomic Mating." Frontiers in Genetics (2016).

Lehermeier at al. "Genetic gain increases by applying the usefulness criterion with improved variance prediction in selection of crosses" Genetics (2017).

Broman et al. "R/qtl: QTL mapping in experimental crosses." Bioinformatics (2003).

VanRaden, Paul M. ''Efficient methods to compute genomic predictions.'' Journal of dairy science (2008).}
\examples{
	\dontrun{
library(GenomicMating)
N=10
nmarkers=200
Markers<-c()
for (i in 1:N){
  Markers<-rbind(Markers,rbinom(nmarkers, 2,.1)-1)
}

Markers2<-c()
for (i in 1:N){
  Markers2<-rbind(Markers2,rbinom(nmarkers, 2,.1)-1)
}

markereffects<-rep(0,nmarkers)
markereffects[sample(1:nmarkers,nmarkers/2)]<-rnorm(nmarkers/2)
Markers[1:5,1:5]
library(parallel)
K=Amat.pieces(rbind(Markers), pieces=5) 

K2=Amat.pieces(rbind(Markers,Markers2), pieces=5) 
K[1:5,1:5]
rownames(Markers)<-paste("l", 1:nrow(Markers),sep="_")
rownames(Markers2)<-paste("l", (nrow(Markers)+1):(nrow(Markers)+nrow(Markers2)),sep="_")
rownames(K2)<-colnames(K2)<-c(rownames(Markers),rownames(Markers2))
rownames(K)<-colnames(K)<-c(rownames(Markers))

markereffects<-rep(0,nmarkers)
markereffects[sample(1:nmarkers,nmarkers/2)]<-rnorm(nmarkers/2)
markereffects2<-rep(0,nmarkers)
markereffects2[sample(1:nmarkers,nmarkers/2)]<-rnorm(nmarkers/2)
markereffects3<-rep(0,nmarkers)
markereffects3[sample(1:nmarkers,nmarkers/2)]<-rnorm(nmarkers/2)

markermap=as.matrix(data.frame(chr=rep(1,nmarkers),pos=seq(0,1,length=nmarkers)))

map<-cbind(1:nmarkers,1,seq(0,1e+2, length=nmarkers))
map<-qtl::sim.map(len=c(.5), n.mar=nmarkers, anchor.tel=TRUE,
  include.x=FALSE, sex.sp=FALSE, eq.spacing=FALSE)
map<-cbind(1:nmarkers,1,map[[1]])
dim(map)

rownames(K)<-colnames(K)<-rownames(Markers)<-1:nrow(Markers)
rownames(map)<-1:ncol(Markers)
sum(is.na(map))

gasols5<-getGaSolutionsFrontierMultiTraitSimcross(Markers=Markers, 
K=K,map=map, markereffectslist=list(markereffects, markereffects2),
nmates=10,npopGA=10, nitGA=10,mc.cores=1,mutprob=0.999, 
nSim = 10,simtype="riself")

gasols5[[1]]
pairs(gasols5[[1]])
}
}
\author{Deniz Akdemir, Julio Isidro Sanch\'ez, Hanna Haikka, Itaraju Baracuhy Brum}