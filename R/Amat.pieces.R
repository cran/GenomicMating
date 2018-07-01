Amat.pieces<-function(M, pieces=10, mc.cores=1){
  #M is coded as -1,0,1; no missing
  
AmatPieces<-function(M){
  pvec<-matrix(apply(M, 2, function(x){mean(x+1)/2}), ncol=1)
  MMt<-tcrossprod(M+1-2*matrix(1, nrow=nrow(M), ncol=1)%*%t(pvec))
  return(list(pvec, MMt))
}


CombAmatPieces<-function(Amatpiecesout){
  nparts<-length(Amatpiecesout)
  NumAmat<-Amatpiecesout[[1]][[2]]
  denomAmat<-sum(2*Amatpiecesout[[1]][[1]]*(1-Amatpiecesout[[1]][[1]]))
  for (i in 2:length(Amatpiecesout)){
    NumAmat<-NumAmat+Amatpiecesout[[i]][[2]]
    denomAmat<-denomAmat+sum(2*Amatpiecesout[[i]][[1]]*(1-Amatpiecesout[[i]][[1]]))
    
  }
  return(NumAmat/denomAmat)
}

x <- 1:ncol(M)
n <- pieces
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

listforlapply<-chunk(x,n)

lapplyfunc<-function(x){
  return(AmatPieces(M[,x]))
}
lapplyout<-mclapply(X=listforlapply, FUN=lapplyfunc, mc.cores=mc.cores)


A<-CombAmatPieces(lapplyout)

return(A)

}
 


