\name{getGaSolutions}
\alias{getGaSolutions}
\alias{Kmatfunc}
\alias{calculatecrossvalue} 
\alias{getstatsM1}
\alias{getstatsM2}
\alias{getstatsM3}
\alias{getstatsfromsim}
\alias{pairs3d}
\alias{par.name}
\alias{mapfunct}
\alias{par.position}
\alias{tails}

\title{getGaSolutions}
\description{Genomic mating is similar to genomic selection in terms of estimating marker effects, but in genomic mating the genetic information and the estimated marker effects are used to decide which genotypes should be crossed to obtain the next breeding population. This program uses genetic algorithm to obtain a solution that minimizes inbreeding and maximize gain and usefulness for a given set of importance weights. 
}
\usage{getGaSolutions(Markers, Markers2=NULL,K, markereffects,
markermap=NULL, nmates=NULL,minparents=10, impinbreedstepsize=.02,
impvar=.1, impforinbreed=.7,keepbest=TRUE, npopGA=100, nitGA=100,
plotiters=TRUE,nelite=10, mutprob=1, mc.cores=1, miniters=100,
minitbefstop=80, tolparconv=1e-6, noself=F,method=1L, type=0,
generation=1,plotMates=TRUE)
}
\arguments{
\item{Markers}{The matrix of markers rows corresponding to individuals and columns for markers, the markers scores are coded as -1,0,1. (For Method=3 the markers are coded as probabilities between 0 and 1.)}
 \item{Markers2}{The matrix of markers rows corresponding to individuals and columns for markers, the markers scores are coded as -1,0,1. (For Method=3 the markers are coded as probabilities between 0 and 1.)}
\item{K}{symmetric genomic relationship matrix, the order of the row and columns of this matrix should follow the order of genotypes in the rows of \code{rbind(Markers, Markers2)}.}
 \item{markermap}{A map for markers. two columns, firts column is named chr second named pos for the chromosome and position of the markers specified above}
\item{markereffects}{effects of markers for a trait}
\item{nmates}{number of mates to select, default value is NULL (number of mates is equal to number of mates)}
\item{minparents}{minimum number of parents in the solution (importance parameter for inbreeding is increased till minimum number of parents are included in the mating solution), minimum is 1.}
\item{impinbreedstepsize}{stepsize for importance parameter for inbreeding to be increased till minimum number of parents are included in the mating solution,}
\item{impvar}{importance parameter of the cross variance term}

 \item{impforinbreed}{importance parameter for inbreeding}
 \item{keepbest}{logical. Default value is TRUE. GA parameter, wheteher to keep the best solution in each iteration.}
\item{npopGA}{genetic algorithm parameter: number of solutions generated at each cycle of the GA}
\item{nitGA}{genetic algorithm parameter: number of GA cycles before algorithm stops }
\item{plotiters}{genetic algorithm parameter: if TRUE the value of the objective function over iterations will be plotted}
\item{nelite}{genetic algorithm parameter:number of elite solutions selected at each cycle of the GA}
\item{mutprob}{genetic algorithm parameter: mutation probability}
\item{mc.cores}{genetic algorithm parameter: number of cores to use}
\item{miniters}{genetic algorithm parameter: minimum number of GA cycles before algorithm stops }
\item{minitbefstop}{genetic algorithm parameter: minimum number of GA cycles before algorithm continues when the tolerance is reached (no change in the criterion value)}
\item{tolparconv}{genetic algorithm parameter: the maximum change in criterion value accepted for conbergence.}
\item{noself}{Is selfing allowed? (TRUE or FALSE)}
\item{method}{Which method to use? (1,2,3) See Detalis.}
\item{type}{Only for method=2. Type of crosses (1 (DH), 2 (RISELF)).}
\item{generation}{Only for method=2. Generation at which the cross variances are calculated.}
 \item{plotMates}{Plot a final 3D plot of solutions. (TRUE or FALSE)}
}
\value{Returns a list with three elements: the first element in this list is the list of mates in the best solution, the second element in the list is the criterion values for the best solutions through the iterations. The last item in the list is a list itself that contains the values of the statistics Gain, Usefulness and Inbreeding.}
\details{The efficient mating problem can be stated as an optimization problem as follows:
minimize w.r.t. \eqn{P_{32}} \eqn{ r( \lambda_1, \lambda_2, P_{32}) = -(1-\lambda_1-\lambda_2)*Gain(P_{32}) -\lambda_1*Usefulness(P_{32}) + \lambda_2 * Inbreeding(P_{32})} where \eqn{0\leq\lambda_1, \lambda_2\leq 1} ,\eqn{0\leq\lambda_2, \lambda_2\leq 1,} \eqn{\lambda_1+\lambda_2\leq 1} and the minimization is over the space of the mating matrices \eqn{P_{32}} construction of which is described in detail in the cited article.

Gain(P) for a mating design P is calculated as $P g$ where g is the vector of genomically estimated breeding values for the parents.
Inbreeding(P) for a mating design P is calculated as $trace(P K P'+D(P))$ where K is the matrix of genomically estimated relationship matrix for the parents and D(P) is a diagonal matrix for adjustment of the parent relationship matrix to progeny relationship matrix. 

Usefulness(P) measures the variance of a mate pair. The average or the sum of usefulnesses for all pairs in a mating plan can be used to measure the usefulness of a mating plan. 
There are three options for the calculation of usefulness. Method 1 uses the calculations in ''Efficient Breeding by Genomic Mating'', Method 2 uses the calculations in ''Genetic gain increases by applying the usefulness criterion with improved variance prediction in selection of crosses'' without the estimation variance terms. Method 2 comes with two types (DH (type=0) or riself (type=1)) and each of these types can be applied for progeny at a specified ''generation''. Method 3 is for polyploid organisms, where the marker data is recorded as proportions ofalleles at genomewide loci. 


}
\references{
Akdemir & Sanchez. "Efficient Breeding by Genomic Mating." Frontiers in Genetics (2016).

Lehermeier at al. "Genetic gain increases by applying the usefulness criterion with improved variance prediction in selection of crosses" Genetics (2017).

Broman et al. "R/qtl: QTL mapping in experimental crosses." Bioinformatics (2003).

VanRaden, Paul M. ''Efficient methods to compute genomic predictions.'' Journal of dairy science (2008).
}
\examples{
\dontrun{
library(GenomicMating)

###Create 100 markers for two sets of populations of size 20.
N=20
nmarkers=100
Markers<-c()
for (i in 1:N){
  Markers<-rbind(Markers,rbinom(nmarkers, 2,.1)-1)
}


Markers2<-c()
for (i in 1:N){
  Markers2<-rbind(Markers2,rbinom(nmarkers, 2,.1)-1)
}

###Marker effects for a trait.
markereffects<-rep(0,nmarkers)
markereffects[sample(1:nmarkers,nmarkers/2)]<-rnorm(nmarkers/2)
Markers[1:5,1:5]

#######Relationship matrices (K only for the first population.
##K2 for both populations together.)
#library(parallel)
K=Amat.pieces(rbind(Markers), pieces=5) 

K2=Amat.pieces(rbind(Markers,Markers2), pieces=5) 
K[1:5,1:5]

####putting names
rownames(Markers)<-paste("l", 1:nrow(Markers),sep="_")
rownames(Markers2)<-paste("l", (nrow(Markers)+1):(nrow(Markers)+
nrow(Markers2)),sep="_")
rownames(K2)<-colnames(K2)<-c(rownames(Markers),rownames(Markers2))
rownames(K)<-colnames(K)<-c(rownames(Markers))


###Best genotype in pop 1
which.max(Markers\%*\%markereffects)
markermap=as.matrix(data.frame(chr=rep(1,nmarkers),
pos=seq(0,1,length=nmarkers)))

colnames(Markers)<-1:nmarkers

########Mating within pop 1, using method 1. 
########Adjust genetic algorithm paparmeters for convergence.

gasols<-getGaSolutions(Markers=Markers,Markers2=NULL, K=K,
markereffects=markereffects,markermap=markermap,nmates=10,
minparents=3, impinbreedstepsize=.02, impvar=.01, 
impforinbreed=.01,npopGA=50, nitGA=10, miniters=10,minitbefstop=20,
plotiters=TRUE,mc.cores=1,nelite=20, mutprob=0.8, noself=TRUE,
method=1, type=0L, generation=0L)

gasols


######Mating between pop1 and pop2. Method 1.

gasols1<-getGaSolutions(Markers=Markers,Markers2=Markers2, K=K2,
markereffects,markermap=markermap,nmates=10,
    minparents=3, 
    impinbreedstepsize=.02, impvar=.02, 
    impforinbreed=.07,
    npopGA=50, nitGA=10, miniters=10,minitbefstop=20,
    plotiters=TRUE,
    mc.cores=2,nelite=20, mutprob=0.8, noself=F, method=1,
    type=0L, generation=0L)


######Mating between pop1 and pop2. Method 2.

gasols2<-getGaSolutions(Markers=Markers,Markers2=Markers2, K=K2,
markereffects,markermap=markermap,nmates=10,
    minparents=3, 
    impinbreedstepsize=.02, impvar=.02, 
    impforinbreed=.07,
    npopGA=50, nitGA=10, miniters=10,minitbefstop=20,
    plotiters=TRUE,
    mc.cores=2,nelite=20, mutprob=0.8, noself=F, method=2,
    type=0L, generation=0L)

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
rownames(Markers2)<-paste("l", (nrow(Markers)+1):(nrow(Markers)+nrow(Markers2)),sep="_")
rownames(K2)<-colnames(K2)<-c(rownames(Markers),rownames(Markers2))
rownames(K)<-colnames(K)<-c(rownames(Markers))

which.max(Markers\%*\%markereffects)
markermap=as.matrix(data.frame(chr=rep(1,nmarkers),pos=seq(0,1,length=nmarkers)))

colnames(Markers)<-1:nmarkers


gasols3<-getGaSolutions(Markers=Markers,Markers2=Markers2, K=K2,
markereffects,markermap=markermap,nmates=10,
    minparents=1, 
    impinbreedstepsize=.02, impvar=.02, 
    impforinbreed=.07,
    npopGA=50, nitGA=10, miniters=10,minitbefstop=20,plotiters=TRUE,
    mc.cores=1,nelite=20, mutprob=0.8, noself=F, method=3,
    type=0L, generation=0L)


gasols3

}

}
\author{Deniz Akdemir, Julio Isidro Sanch\'ez, Hanna Haikka, Itaraju Baracuhy Brum}
