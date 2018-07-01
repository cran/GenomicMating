\name{plotGM}
\alias{plotGM}

\title{plotGM}
\description{For plotting GA results. See examples.
}
\usage{plotGM(GMsols, type="3D", plotly = FALSE, idealsol=NULL, traitnum=1)
}
\arguments{
\item{GMsols}{Output of getGaSolutionsMultiTrait of getGaSolutionsMultiTraitSinCross.}
 \item{type}{Options are "3D", "SOM", "SOM2".}
 \item{plotly}{Logical, default is FALSE. Uses plotly for 3D plot if plotly is inslalled.}
\item{idealsol}{For coloring the plot. Defaults is NULL. Otherwise, a vector of the same length as the GMsols statistics.}
\item{traitnum}{which trait (order of trait in the markereffectslist).}
}
\value{NULL}
\details{See examples}
\references{Akdemir, Deniz, and Julio I. Sanchez. ''Efficient Breeding by Genomic Mating.'' Frontiers in Genetics 7 (2016).}
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
rownames(Markers2)<-paste("l", (nrow(Markers)+1):(nrow(Markers)+nrow(Markers2)),sep="_")
rownames(K2)<-colnames(K2)<-c(rownames(Markers),rownames(Markers2))
rownames(K)<-colnames(K)<-c(rownames(Markers))

which.max(Markers\%*\%markereffects)
markermap=as.matrix(data.frame(chr=rep(1,nmarkers),pos=seq(0,1,length=nmarkers)))

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

