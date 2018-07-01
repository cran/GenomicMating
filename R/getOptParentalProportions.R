getOptParentalProportions<-function(Amat, gebvs, lambda, ul){
  N<-nrow(Amat)
  H<-2*lambda*Amat
  d<--(1-lambda)*gebvs
  A<-matrix(1,nrow=1,ncol=N)
  u<-matrix(1, ncol=1,nrow=N)*ul
  log <- capture.output({
     solutionqp<-suppressMessages(invisible(LowRankQP::LowRankQP(Vmat=H,dvec=d,Amat=A,bvec=1,uvec=u,method="LU", verbose=F)))
   })
  out<-c(solutionqp$alpha, lambda,crossprod(solutionqp$alpha,gebvs),crossprod(solutionqp$alpha,Amat%*%solutionqp$alpha),crossprod(solutionqp$alpha,gebvs)/crossprod(solutionqp$alpha,Amat%*%solutionqp$alpha))
  if (is.null(rownames(Amat))){namesAmat<-paste("l",1:N,sep="")}else{namesAmat<-rownames(Amat)}
  names(out)<-c(namesAmat, "lambda", "Gain","Inbreeding","G/I ratio")
  return(out)
}




plotOPFrontier<-function(Amat, gebvs, lambdavec=seq(1e-5,1-1e-5, length=100), ul=1, mc.cores=1, identify=FALSE){
  
  out<-suppressMessages(suppressWarnings(t(simplify2array(mclapply(lambdavec, function(lambda){getOptParentalProportions(Amat, gebvs, lambda, ul)}, mc.cores=mc.cores)))))
  N<-nrow(Amat)
  if (is.null(rownames(Amat))){namesAmat<-paste("l",1:N,sep="")}else{namesAmat<-rownames(Amat)}
  colnames(out)<-c(namesAmat, "lambda", "Gain","Inbreeding","G/I ratio")
  ScaledIandG<-scale(cbind(out[,ncol(out)-2], out[,ncol(out)-1]))
  max1<-max(ScaledIandG[,1])
  min2<-min(ScaledIandG[,2])
  distD<-(as.matrix(dist(rbind(c(max1,min2),ScaledIandG))))[1,-1]
  rbPal <- colorRampPalette(c('red','blue'))
  #This adds a column of color values
  # based on the y values
  Col <- rbPal(10)[distD]
  plot(out[,ncol(out)-2], out[,ncol(out)-1], xlab="Gain", ylab="Inbreeding", col=Col)
  rownames(out)<-paste("sol_", 1:nrow(out), sep="")
  if (identify){identify(out[,ncol(out)-2], out[,ncol(out)-1], labels=row.names(out), plot=TRUE, col="blue")}
  return(out)
}
