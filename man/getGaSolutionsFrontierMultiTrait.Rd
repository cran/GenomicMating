\name{getGaSolutionsFrontierMultiTrait}
\alias{getGaSolutionsFrontierMultiTrait}
\title{getGaSolutionsFrontierMultiTrait}
\description{Generalizes the approach in getGaSolutions to multiple traits specified by a list of marker effects.}
\usage{getGaSolutionsFrontierMultiTrait(Markers,Markers2=NULL, K,
markereffectslist,markermap,nmates=NULL,npopGA, nitGA, mutprob,
mc.cores, noself=F,method=1, type=0L, generation=0L, plotiters=F)}
\arguments{
\item{Markers}{The matrix of markers rows corresponding to individuals and columns for markers, the markers scores are coded as -1,0,1. (For Method=3 the markers are coded as probabilities between 0 and 1.)}
 \item{Markers2}{The matrix of markers rows corresponding to individuals and columns for markers, the markers scores are coded as -1,0,1. (For Method=3 the markers are coded as probabilities between 0 and 1.)}
\item{K}{symmetric genomic relationship matrix, the order of the row and columns of this matrix should follow the order of genotypes in the rows of \code{rbind(Markers, Markers2)}.}
  \item{markermap}{a map for markers. two columns, firts column is named chr second named pos for the chromosome and position of the markers specified above}
 \item{markereffectslist}{effects of markers for several traits given as a list}
  \item{nmates}{number of mates to select, default value is NULL (number of mates is equal to number of genotypes)}
\item{npopGA}{genetic algorithm parameter: number of solutions generated at each cycle of the GA}
\item{nitGA}{genetic algorithm parameter: number of GA cycles before algorithm stops }
\item{mutprob}{genetic algorithm parameter: mutation probability}
\item{mc.cores}{genetic algorithm parameter: number of cores to use}
\item{noself}{Is selfing allowed? (TRUE or FALSE)}
\item{method}{Which method to use? (1,2,3) See Details.}
\item{type}{Only for method=2. Type of crosses (0 (DH), 1 (RISELF)).}
\item{generation}{Only for method=2. Generation at which the cross variances are calculated.}
\item{plotiters}{genetic algorithm parameter: if TRUE the value of the objective function over iterations will be plotted}}
\value{Returns a list with two elements: the first element in this list is a list of solutions found on the frontier, the second element is the matrix of criterion values (Gain, Usefulness, and Inbreeding) corresponding to these solutions.}
\details{This program uses genetic algorithm to produce a number of solutions on the frontier curve simultaneously for the multi-objective optimization problem which is defined by minimization of \eqn{-Gain_j(P_{32}),} \eqn{-Usefulness_j(P_{32})} for \eqn{j=1,2,...,ntraits} and \eqn{Inbreeding(P_{32})} with respect to \eqn{P_{32}.}

Gain(P) for a mating design P is calculated as \eqn{P g} where g is the vector of genomically estimated breeding values for the parents.

Inbreeding(P) for a mating design P is calculated as \eqn{trace(P K P'+D(P))} where K is the matrix of genomically estimated relationship matrix for the parents and D(P) is a diagonal matrix for adjustment of the parent relationship matrix to progeny relationship matrix. 

Usefulness(P) measures the variance of a mate pair. The average or the sum of usefulnesses for all pairs in a mating plan can be used to measure the usefulness of a mating plan. 
There are three options for the calculation of usefulness. Method 1 uses the calculations in ''Efficient Breeding by Genomic Mating'', Method 2 uses the calculations in ''Genetic gain increases by applying the usefulness criterion with improved variance prediction in selection of crosses'' without the estimation variance terms. Method 2 comes with two types (DH (type=0) or riself (type=1)) and each of these types can be applied for progeny at a specified ''generation''. Method 3 is for polyploid organisms, where the marker data is recorded as proportions ofalleles at genomewide loci.}
\references{
Akdemir & Sanchez. "Efficient Breeding by Genomic Mating." Frontiers in Genetics (2016).

Lehermeier at al. "Genetic gain increases by applying the usefulness criterion with improved variance prediction in selection of crosses" Genetics (2017).

Broman et al. "R/qtl: QTL mapping in experimental crosses." Bioinformatics (2003).

VanRaden, Paul M. ''Efficient methods to compute genomic predictions.'' Journal of dairy science (2008).
}

\examples{
	\dontrun{
library("GenomicMating")

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

K=Amat.pieces(rbind(Markers), pieces=5) 
K2=Amat.pieces(rbind(Markers,Markers2), pieces=5) 

rownames(Markers)<-paste("l", 1:nrow(Markers),sep="_")
rownames(Markers2)<-paste("l", (nrow(Markers)+1):(nrow(Markers)+nrow(Markers2)),sep="_")

rownames(K2)<-colnames(K2)<-c(rownames(Markers),rownames(Markers2))
rownames(K)<-colnames(K)<-c(rownames(Markers))

#######Two sets of marker effects

markereffects<-rep(0,nmarkers)
markereffects[sample(1:nmarkers,nmarkers/2)]<-rnorm(nmarkers/2)
markereffects2<-rep(0,nmarkers)
markereffects2[sample(1:nmarkers,nmarkers/2)]<-rnorm(nmarkers/2)


gasols4<-getGaSolutionsFrontierMultiTrait(Markers=Markers,
Markers2=Markers2,K=K2,
markereffectslist=list(markereffects,markereffects2),
markermap=markermap,nmates=20,npopGA=100, nitGA=10, 
mc.cores=1, mutprob=0.99,method=2,
type=0, generation=3, plotiters= TRUE)

str(gasols4)
gasols4[[1]][1:5,]
}
}
\author{Deniz Akdemir, Julio Isidro Sanch\'ez, Hanna Haikka, Itaraju Baracuhy Brum}