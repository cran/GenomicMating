getGaSolutions<-function(Markers, K, markereffects,nmates=NULL,minparents=10, impinbreedstepsize=.02, impvar=.1, impforinbreed=.7,npopGA=100, nitGA=100, plotiters=TRUE,nelite=10, mutprob=1, mc.cores=1, miniters=100,minitbefstop=80, tolparconv=1e-6){

 
  if (minitbefstop>miniters){minitbefstop<-miniters-1}
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
  
  getstatsout<-getstats(Markers, K, markereffects, P=Z)
 
  return(-(1-impvar-impforinbreed)*getstatsout[1]-impvar*getstatsout[3]+impforinbreed*getstatsout[2])
}

GAobjfunclist<-function(indexesinallcombs){
  Z<-matrix(0,nrow=nmates, ncol=N)
  for (i in 1:length(indexesinallcombs)){
    Z[i,allcombs[indexesinallcombs[i],]]<-1/2
    if (allcombs[indexesinallcombs[i],][1]==allcombs[indexesinallcombs[i],][2]){Z[i,allcombs[indexesinallcombs[i],]]<-1}
  }
  
  getstatsout<-getstats(Markers, K, markereffects, P=Z)
  
  return(list(getstatsout[1],getstatsout[3],getstatsout[2]))
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
    if (iters>miniters){
      if (length(table(round(meanvec[((iters-minitbefstop):iters)],tolparconv)))==1){
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

nameslines<-rownames(Markers)

if (!is.null(nameslines)){ if (length(nameslines)==nrow(Markers)){
  solGA<-cbind(nameslines[solGA[,1]],nameslines[solGA[,2]])
}}
stats<-GAobjfunclist(solutionGA)
return(list(solGA, stat,stats))
}










######################################################
getGaSolutionsFrontier<-function(Markers, K, markereffects,nmates=NULL,npopGA, nitGA, mutprob, mc.cores){
  ##############################list of parents optimization
  N=nrow(Markers)
  if (is.null(nmates)){nmates=N}
  allcombs<-matrix(nrow=0,ncol=2)
  ij=1
  for (i in 1:N){for(j in i:N){allcombs<-rbind(allcombs,c(i,j));ij+ij+1}}
  
    GAobjfunc<-function(indexesinallcombs){
      Z<-matrix(0,nrow=nmates, ncol=N)
      for (i in 1:length(indexesinallcombs)){
        Z[i,allcombs[indexesinallcombs[i],]]<-1/2
        if (allcombs[indexesinallcombs[i],][1]==allcombs[indexesinallcombs[i],][2]){Z[i,allcombs[indexesinallcombs[i],]]<-1}
      }
      
      getstatsout<-getstats(Markers, K, markereffects, P=Z)
      out<-list(-getstatsout[1],-getstatsout[3],getstatsout[2])
      return(out)
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
    
    
    GAfunc<-function (ntoselect, npop,  mutprob, 
                      niterations) 
    {
      
      Candidates<-1:nrow(allcombs)
      
      InitPop <- mclapply(1:npop, function(x) {
        return(sample(Candidates, ntoselect, replace=T))
      }, mc.cores=mc.cores)
      
      InitPopFuncValues <- matrix(as.numeric(unlist(mclapply(InitPop,FUN = function(x) {
        GAobjfunc(x)}, mc.cores = mc.cores, mc.preschedule = T))), ncol= 3, byrow=T)
      
      
      frontier3<- which(!is_dominated(t(InitPopFuncValues[,1:3])))
      xy.f <- InitPopFuncValues[frontier3, ]
      scatterplot3d(xy.f[,c(1,2,3)], highlight.3d=T)
      
      
      #
      # Visualization.
      #
     # plot(xy.f[,c(1,2)], xlab="X", ylab="Y", pch=19, main="Quasiconvex Hull", col="red")
     # points(InitPopFuncValues[,c(1,2)], col="blue")
      
      ElitePop <- mclapply(frontier3, FUN = function(x) {
        return(InitPop[[x]])
      }, mc.cores = mc.cores, mc.preschedule = T)
      ElitePopFuncValues <- InitPopFuncValues[frontier3,]
      
      for (iters in 1:niterations) {
        
        CurrentPop <- GenerateCrossesfromElites(Elites = ElitePop, 
                                                Candidates = Candidates, npop = npop, mutprob = mutprob)
        CurrentPop<-c(CurrentPop, ElitePop[1])
        CurrentPopFuncValues <- matrix(as.numeric(unlist(mclapply(CurrentPop,
                                                           FUN = function(x) {
                                                             GAobjfunc(x)
                                                           }, mc.cores=mc.cores, mc.preschedule = T))), ncol= 3, byrow=T)
        
        
        CurrentPop <- c(CurrentPop,ElitePop)
        notduplicated<-!duplicated(CurrentPop)
        CurrentPop<-CurrentPop[notduplicated]
        CurrentPopFuncValues<-rbind(CurrentPopFuncValues,ElitePopFuncValues)
        CurrentPopFuncValues<-CurrentPopFuncValues[notduplicated,]
        # frontier <- hull.qc(CurrentPopFuncValues)
        
        # frontier2<- nondominated_points(t(CurrentPopFuncValues[,1:2]))
        # print(frontier2)
        # str(frontier2)
        frontier3<- which(!is_dominated(t(CurrentPopFuncValues[,1:3])))
        
        xy.f <- CurrentPopFuncValues[frontier3, ]
        
        
        #
        # Visualization.
        #
        scatterplot3d(xy.f[,c(1,2,3)], highlight.3d=T)
        #points3d(CurrentPopFuncValues[,c(1,2,3)], col="blue")
        
        ElitePop <- mclapply(frontier3, FUN = function(x) {
          return(CurrentPop[[x]])
        }, mc.cores = mc.cores, mc.preschedule = T)
        ElitePopFuncValues <- CurrentPopFuncValues[frontier3,]
      }
   
      return(list(ElitePop, ElitePopFuncValues))
    }
    
    outGA<-GAfunc(ntoselect=nmates, npop=npopGA, mutprob=mutprob, 
                  niterations=nitGA) 
    
    
    listofsols<-lapply(outGA[[1]], FUN=function(x){allcombs[x,]})
    nameslines<-rownames(Markers)
    
    if (!is.null(nameslines)){ if (length(nameslines)==nrow(Markers)){
    listofsols<-lapply(listofsols, FUN=function(x){
      return(cbind(nameslines[x[,1]],nameslines[x[,2]]))
    })
    }}
    outGA<-outGA[[2]]
    
    return(list(outGA, listofsols))
    
}
    








