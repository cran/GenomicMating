getGaSolutions<-function(Markers, K, markereffects,nmates=NULL,minparents=10, impinbreedstepsize=.2, impvar, impforinbreed,npopGA, nitGA, plotiters,nelite, mutprob, mc.cores, miniters=100,minitbefstop=80, tolparconv=1e-6){

 converged=FALSE
  ##############################list of parents optimization
N=nrow(Markers)
if (is.null(nmates)){nmates=N}
  allcombs<-matrix(nrow=0,ncol=2)
ij=1
for (i in 1:N){for(j in i:N){allcombs<-rbind(allcombs,c(i,j));ij+ij+1}}

solnnotfound=TRUE
while (((solnnotfound)&&(!converged))){

GAobjfunc<-function(indexesinallcombs){
  Z<-matrix(0,nrow=nmates, ncol=N)
  for (i in 1:length(indexesinallcombs)){
    Z[i,allcombs[indexesinallcombs[i],]]<-1/2
    if (allcombs[indexesinallcombs[i],][1]==allcombs[indexesinallcombs[i],][2]){Z[i,allcombs[indexesinallcombs[i],]]<-1}
  }
  
  getstatsout<-getstats(Markers, K, markereffects, P=Z, impvar=impvar)
 
  return(-getstatsout[3]+impforinbreed*getstatsout[2])
}


makeonecross<-function (x1, x2, Candidates, mutprob = 0.5) 
{
  n1 <- length(unlist(x1))
  n2 <- length(unlist(x2))
  n <- min(c(min(n1), min(n2)))
  x1x2 <- c(unlist(x1), unlist(x2))
  
  cross <- sort(sample(x1x2, n, replace = T))
  randnum <- runif(1)
  if (randnum < mutprob) {
    nmutate=sample(1:2,1)
    cross[sample(1:n, nmutate)] <- sample(setdiff(Candidates, cross), nmutate)
  }
  return(sort(cross))
}

GenerateCrossesfromElites<-function (Elites, Candidates, npop, mutprob) 
{
  newcrosses <- mclapply(1:npop, FUN = function(x) {
    x1 <- Elites[[sample(1:length(Elites), 1)]]
    x2 <- Elites[[sample(1:length(Elites), 1)]]
    return(makeonecross(x1 = x1, x2 = x2, Candidates = Candidates, 
                        mutprob = mutprob))
  }, mc.cores=mc.cores)
  return(newcrosses)
}


GAfunc<-function (ntoselect, npop, nelite, mutprob, 
                  niterations, plotiters = FALSE) 
{
  
  Candidates<-1:nrow(allcombs)
 
  InitPop <- mclapply(1:npop, function(x) {
    return(sample(Candidates, ntoselect, replace=T))
  }, mc.cores=mc.cores)
  
  InitPopFuncValues <- as.numeric(unlist(mclapply(InitPop, FUN = function(x) {
    GAobjfunc(x)
  }, mc.cores=mc.cores)))
  orderofInitPop <- order(InitPopFuncValues, decreasing = FALSE)
  ElitePop <- mclapply(orderofInitPop[1:nelite], FUN = function(x) {
    return(InitPop[[x]])
  }, mc.cores=1)
  ElitePopFuncValues <- InitPopFuncValues[orderofInitPop[1:nelite]]
  meanvec <- c()
  for (iters in 1:niterations) {
   
    CurrentPop <- GenerateCrossesfromElites(Elites = ElitePop, 
                                            Candidates = Candidates, npop = npop, mutprob = mutprob)
    CurrentPop<-c(CurrentPop, ElitePop[1])
    CurrentPopFuncValues <- as.numeric(unlist(mclapply(CurrentPop, 
                                                     FUN = function(x) {
                                                       GAobjfunc(x)
                                                     }, mc.cores=mc.cores)))
   # orderofCurrentPop <- sort(rank(CurrentPopFuncValues, ties.method="min"), decreasing=F)
   orderofCurrentPop <- order((CurrentPopFuncValues), decreasing = FALSE)
    ElitePop <- mclapply(orderofCurrentPop[1:nelite], FUN = function(x) {
      return(CurrentPop[[x]])
    }, mc.cores=mc.cores)
    
    
    ElitePopFuncValues <- CurrentPopFuncValues[orderofCurrentPop[1:nelite]]
    
    meanvec <- c(meanvec, min(ElitePopFuncValues))
    if (plotiters) {
      plot(-meanvec)
    }
    if (iters >miniters){
      if (length(table(round(meanvec[((length(meanvec)-minitbefstop):length(meanvec))],tolparconv)))==1){
       converged=TRUE
         break
      }
    }
  }
  ElitePop[[nelite+1]]<-meanvec
  return(ElitePop)
}

outGA<-GAfunc(ntoselect=nmates, npop=npopGA, nelite=nelite, mutprob=mutprob, 
              niterations=nitGA, plotiters = plotiters) 

solutionGA<-outGA[[1]]
stat<-outGA[[nelite+1]]
solGA<-allcombs[solutionGA,]
if (length(table(unlist(c(solGA))))>=minparents){
  solnnotfound=FALSE
} else {impforinbreed= impforinbreed+impinbreedstepsize}

}

return(list(solGA, stat))
}

