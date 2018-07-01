
############################################

getGaSolutions<-function(Markers, Markers2=NULL,K, markereffects,markermap=NULL, nmates=NULL,minparents=10, impinbreedstepsize=.02, impvar=.1, impforinbreed=.7,keepbest=TRUE, npopGA=100, nitGA=100, plotiters=TRUE,nelite=10, mutprob=1, mc.cores=1, miniters=100,minitbefstop=80, tolparconv=1e-6, noself=F,method=1L, type=0, generation=1,plotMates=TRUE){
  if (minitbefstop>miniters){minitbefstop<-miniters-1}
  converged=FALSE
  ##############################list of parents optimization
  
  if (!is.null(Markers2)){
    N1<-nrow(Markers)
    N2<-nrow(Markers2)
    
    MARKERS<-rbind(Markers,Markers2)
  } else {
    MARKERS=Markers
  }
  
  N=nrow(MARKERS)
  
  if (is.null(Markers2)){
    if (is.null(nmates)){nmates=N}
    allcombs<-matrix(nrow=0,ncol=2)
    ij=1
    if (noself){ 
      for (i in 1:(N-1)){for(j in (i+1):N){allcombs<-rbind(allcombs,c(i,j));ij+ij+1}}
    } else {
      for (i in 1:N){for(j in i:N){allcombs<-rbind(allcombs,c(i,j));ij+ij+1}}
    }
    
  } else {
    if (is.null(nmates)){nmates=N1}
    allcombs<-matrix(nrow=0,ncol=2)
    ij=1
    for (i in 1:N1){for(j in (N1+1):(N1+N2)){allcombs<-rbind(allcombs,c(i,j));ij+ij+1}}
  }
  
  solnnotfound=TRUE
  while (((solnnotfound)&&(!converged))){
    
    GAobjfunc<-function(indexesinallcombs){
      
      Z<-matrix(0,nrow=nmates, ncol=N)
      for (i in 1:length(indexesinallcombs)){
        Z[i,allcombs[indexesinallcombs[i],]]<-1/2
        if (allcombs[indexesinallcombs[i],][1]==allcombs[indexesinallcombs[i],][2]){Z[i,allcombs[indexesinallcombs[i],]]<-1}
      }
      
      if (method==1){
        #arma::vec getstatsM1(arma::mat Markers, arma::mat K, arma::vec markereffects, arma::mat P) {
        
        getstatsout<-getstatsM1(MARKERS, K, markereffects, P=Z)
      } else if (method==2){
        #
        getstatsout<-getstatsM2(MARKERS, K, markereffects, P=Z,markermap, type=type, generation=generation)
      } else if (method==3){
        getstatsout<-getstatsM3(MARKERS, K, markereffects, P=Z)
      } else {
        getstatsout<-getstatsM1(MARKERS, K, markereffects, P=Z)
        
      }      
      return(-(1-impvar-impforinbreed)*getstatsout[2]-impvar*getstatsout[3]+impforinbreed*getstatsout[1])
    }
    
    GAobjfunclist<-function(indexesinallcombs){
      Z<-matrix(0,nrow=nmates, ncol=N)
      for (i in 1:length(indexesinallcombs)){
        Z[i,allcombs[indexesinallcombs[i],]]<-1/2
        if (allcombs[indexesinallcombs[i],][1]==allcombs[indexesinallcombs[i],][2]){Z[i,allcombs[indexesinallcombs[i],]]<-1}
      }
      
      if (method==1){
        #arma::vec getstatsM1(arma::mat Markers, arma::mat K, arma::vec markereffects, arma::mat P) {

        getstatsout<-getstatsM1(MARKERS, K, markereffects, P=Z)
      } else if (method==2){
        #
        getstatsout<-getstatsM2(MARKERS, K, markereffects, P=Z,markermap, type=type, generation=generation)
      } else if (method==3){
        getstatsout<-getstatsM3(MARKERS, K, markereffects, P=Z)
      } else {
        getstatsout<-getstatsM1(MARKERS, K, markereffects, P=Z)
        
      }
      
      
      return(list(getstatsout[1],getstatsout[2],getstatsout[3]))
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
        return(sample(Candidates, ntoselect, replace=TRUE))
      }, mc.cores=mc.cores)
      
      InitPopFuncValues <- as.numeric(unlist(mclapply(InitPop, FUN = function(x) {
        GAobjfunc(x)
      }, mc.cores=mc.cores)))
      orderofInitPop <- order(InitPopFuncValues, decreasing = FALSE)
      ElitePop <- mclapply(orderofInitPop[1:nelite], FUN = function(x) {
        return(InitPop[[x]])
      }, mc.cores=mc.cores)
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
          plot(-meanvec, xlab="iteration", ylab="Statistic")
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
  
  nameslines<-rownames(MARKERS)
  
  if (!is.null(nameslines)){ if (length(nameslines)==nrow(MARKERS)){
    solGA<-cbind(nameslines[solGA[,1]],nameslines[solGA[,2]])
  }}
  stats<-GAobjfunclist(solutionGA)
  names(stats)<-c("I","G","U")
  names(stat)<-paste("StatVal_Iter", 1:length(stat), sep="")
  outlist<-list(solGA, stat,stats)
  names(outlist)<-c("Mates","StatVal_Iter", "Stats")
  
  ########################################
  if(plotMates){
  ebvs=MARKERS%*%markereffects
  names(ebvs)<-rownames(MARKERS)
  
  
  tableslected<-table(factor(unlist(c(outlist$Mates)), levels=rownames(MARKERS)))
  
  allcombs<-c()
  for (i in 1:N){
    for (j in i:N){
      allcombs<-c(allcombs, paste(rownames(MARKERS)[i],rownames(MARKERS)[j], sep="x"))
    }
  }
  
  tableslected2<-table(factor(paste(outlist$Mates[,1],outlist$Mates[,2], sep="x"), levels=allcombs))
  tableslected2<- sort(tableslected2[tableslected2>0], decreasing=TRUE)
  

  
  PCK<-svd(K)
  PCK<-K%*%PCK$v
  s <- scatterplot3d::scatterplot3d(PCK[,1],PCK[,2],ebvs,cex.symbols=log(tableslected+1.3), 
                     highlight.3d = TRUE,  angle = 120,
                     col.axis = "blue", col.grid = "lightblue", cex.axis = 1.3,
                     cex.lab = 1.1, pch = 16,type = "h", 
                     main = paste("N = ", N, ", ImpVar = ", round(impvar,3), ", ImpInbreed = ",round(impforinbreed,3),"\n I:" ,round(outlist$Stats$I,3)," G:" ,round(outlist$Stats$G,3)," U:" ,round(outlist$Stats$U,3),sep=""),xlab="PC1", ylab="PC2", zlab="BV")
  legend(s$xyz.convert(min(PCK[,1])-.1, max(PCK[,2])+.1, mean(ebvs)-5*sd(ebvs)), yjust=-6,xjust=0,
         legend = paste(as.character(names(tableslected2)),as.character((tableslected2)), sep=" : "), cex = .4)
  
  for (i in 1:nrow(PCK)){
    cordtext<- s$xyz.convert(PCK[i,1]+.05,PCK[i,2]+.05,ebvs[i]+1)
    if (tableslected[i]>0){text(cordtext$x, cordtext$y,names(tableslected)[i])}
  }
  ## now draw a line between points 2 and 3
  
  for (i in 1:nrow(outlist$Mates)){
    p2 <- s$xyz.convert(PCK[rownames(PCK)%in%c(outlist$Mates[i,1]),1],PCK[rownames(PCK)%in%c(outlist$Mates[i,1]),2],ebvs[rownames(PCK)%in%c(outlist$Mates[i,1])])
    p3 <- s$xyz.convert(PCK[rownames(PCK)%in%c(outlist$Mates[i,2]),1],PCK[rownames(PCK)%in%c(outlist$Mates[i,2]),2],ebvs[rownames(PCK)%in%c(outlist$Mates[i,2])])
    segments(p2$x,p2$y,p3$x,p3$y,lwd=2,col=1)
  }
  
  
  
  }
  ##########################################
  return(outlist)
}





######################################################
getGaSolutionsFrontier<-function(Markers, Markers2=NULL,K, markereffects,markermap=NULL,nmates=NULL,npopGA=100, nitGA=100, mutprob=1, mc.cores=1, noself=F,method=1, type=0, generation=1L, plotiters=F){
  ##############################list of parents optimization
  
  if (!is.null(Markers2)){
    N1<-nrow(Markers)
    N2<-nrow(Markers2)
    
    MARKERS<-rbind(Markers,Markers2)
  } else {
    MARKERS=Markers
  }
  
  N=nrow(MARKERS)
  
  if (is.null(Markers2)){
    if (is.null(nmates)){nmates=N}
    allcombs<-matrix(nrow=0,ncol=2)
    ij=1
    if (noself){ 
      for (i in 1:(N-1)){for(j in (i+1):N){allcombs<-rbind(allcombs,c(i,j));ij+ij+1}}
    } else {
      for (i in 1:N){for(j in i:N){allcombs<-rbind(allcombs,c(i,j));ij+ij+1}}
    }
    
  } else {
    if (is.null(nmates)){nmates=N1}
    allcombs<-matrix(nrow=0,ncol=2)
    ij=1
    for (i in 1:N1){for(j in (N1+1):(N1+N2)){allcombs<-rbind(allcombs,c(i,j));ij+ij+1}}
  }
  
  GAobjfunc<-function(indexesinallcombs){
    Z<-matrix(0,nrow=nmates, ncol=N)
    for (i in 1:length(indexesinallcombs)){
      Z[i,allcombs[indexesinallcombs[i],]]<-1/2
      if (allcombs[indexesinallcombs[i],][1]==allcombs[indexesinallcombs[i],][2]){Z[i,allcombs[indexesinallcombs[i],]]<-1}
    }
    
    
    if (method==1){
      #arma::vec getstatsM1(arma::mat Markers, arma::mat K, arma::vec markereffects, arma::mat P) {
      
      getstatsout<-getstatsM1(MARKERS, K, markereffects, P=Z)
    } else if (method==2){
      #
      getstatsout<-getstatsM2(MARKERS, K, markereffects, P=Z,markermap, type=type, generation=generation)
    } else if (method==3){
      getstatsout<-getstatsM3(MARKERS, K, markereffects, P=Z)
    } else {
      getstatsout<-getstatsM1(MARKERS, K, markereffects, P=Z)
      
    }
    
    out<-list(getstatsout[1],-getstatsout[2],-getstatsout[3])
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
                    niterations, type=type, generation=generation) 
  {
    
    Candidates<-1:nrow(allcombs)
    
    InitPop <- mclapply(1:npop, function(x) {
      return(sample(Candidates, ntoselect, replace=TRUE))
    }, mc.cores=mc.cores)
    
    InitPopFuncValues <- matrix(as.numeric(unlist(mclapply(InitPop,FUN = function(x) {
      GAobjfunc(x)}, mc.cores = mc.cores, mc.preschedule = T))), ncol= 3, byrow=TRUE)
    
    
    frontier3<- which(!emoa::is_dominated(t(InitPopFuncValues[,1:3])))
    xy.f <- InitPopFuncValues[frontier3, ]
    if (plotiters){scatterplot3d::scatterplot3d(xy.f[,c(1,2,3)], highlight.3d=TRUE, xlab="I", ylab="-G", zlab="-U")}
    
    
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
                                                                }, mc.cores=mc.cores, mc.preschedule = T))), ncol= 3, byrow=TRUE)
      
      
      CurrentPop <- c(CurrentPop,ElitePop)
      notduplicated<-!duplicated(CurrentPop)
      CurrentPop<-CurrentPop[notduplicated]
      CurrentPopFuncValues<-rbind(CurrentPopFuncValues,ElitePopFuncValues)
      CurrentPopFuncValues<-CurrentPopFuncValues[notduplicated,]
      
      frontier3<- which(!emoa::is_dominated(t(CurrentPopFuncValues[,1:3])))
      xy.f <- CurrentPopFuncValues[frontier3, ]
      
      
      #
      # Visualization.
      #
      if (plotiters){scatterplot3d::scatterplot3d(xy.f[,c(1,2,3)], highlight.3d=TRUE, xlab="I", ylab="-G", zlab="-U")}
      #points3d(CurrentPopFuncValues[,c(1,2,3)], col="blue")
      
      ElitePop <- mclapply(frontier3, FUN = function(x) {
        return(CurrentPop[[x]])
      }, mc.cores = mc.cores, mc.preschedule = T)
      ElitePopFuncValues <- CurrentPopFuncValues[frontier3,]
    }
    
    return(list(ElitePop, ElitePopFuncValues))
  }
  
  outGA<-GAfunc(ntoselect=nmates, npop=npopGA, mutprob=mutprob, 
                niterations=nitGA, type=type, generation=generation) 
  
  
  listofsols<-lapply(outGA[[1]], FUN=function(x){allcombs[x,]})
  nameslines<-rownames(MARKERS)
  
  if (!is.null(nameslines)){ if (length(nameslines)==nrow(MARKERS)){
    listofsols<-lapply(listofsols, FUN=function(x){
      return(cbind(nameslines[x[,1]],nameslines[x[,2]]))
    })
  }}
  outGA<-outGA[[2]]
  
  return(list(outGA, listofsols))
  
}






######################################################
getGaSolutionsFrontierMultiTrait<-function(Markers,Markers2=NULL, K, markereffectslist,markermap,nmates=NULL,npopGA, nitGA, mutprob, mc.cores, noself=F,method=1,
                                           type=0L, generation=0L, plotiters=F){
  ##############################list of parents optimization
  
  
  if(is.null(names(markereffectslist))){
    
    names(markereffectslist)<-paste("T_", 1:length(markereffectslist),sep="")
  }
  colnamesforoutput<-c("I",paste(rep(c("-G_", "-U_"),length(names(markereffectslist))),rep(names(markereffectslist), each=length(markereffectslist)),sep=""))  
  
  if (!is.null(Markers2)){
    N1<-nrow(Markers)
    N2<-nrow(Markers2)
    
    MARKERS<-rbind(Markers,Markers2)
  } else {
    MARKERS=Markers
  }
  
  N=nrow(MARKERS)
  
  if (is.null(Markers2)){
    if (is.null(nmates)){nmates=N}
    allcombs<-matrix(nrow=0,ncol=2)
    ij=1
    if (noself){ 
      for (i in 1:(N-1)){for(j in (i+1):N){allcombs<-rbind(allcombs,c(i,j));ij+ij+1}}
    } else {
      for (i in 1:N){for(j in i:N){allcombs<-rbind(allcombs,c(i,j));ij+ij+1}}
    }
    
  } else {
    if (is.null(nmates)){nmates=N1}
    allcombs<-matrix(nrow=0,ncol=2)
    ij=1
    for (i in 1:N1){for(j in (N1+1):(N1+N2)){allcombs<-rbind(allcombs,c(i,j));ij+ij+1}}
  }
  
  GAobjfunc<-function(indexesinallcombs){
    Z<-matrix(0,nrow=nmates, ncol=N)
    for (i in 1:length(indexesinallcombs)){
      Z[i,allcombs[indexesinallcombs[i],]]<-1/2
      if (allcombs[indexesinallcombs[i],][1]==allcombs[indexesinallcombs[i],][2]){Z[i,allcombs[indexesinallcombs[i],]]<-1}
    }
    
    out<-vector(mode="list")
    
    if (method==1){
      getstatsout<-getstatsM1(MARKERS, K, markereffectslist[[1]], P=Z)
    } else if (method==2){
      getstatsout<-getstatsM2(MARKERS, K, markereffectslist[[1]], P=Z,markermap, type=type, generation=generation)
    } else if (method==3){
      getstatsout<-getstatsM3(MARKERS, K, markereffectslist[[1]], P=Z)
    } else {
      getstatsout<-getstatsM1(MARKERS, K, markereffectslist[[1]], P=Z)
    }
    
    
    
    out<-list(getstatsout[1],-getstatsout[2],-getstatsout[3])
    for (i in 2:length(markereffectslist)){
      if (method==1){
        getstatsout<-getstatsM1(MARKERS, K, markereffectslist[[i]], P=Z)
      } else if (method==2){
     
        getstatsout<-getstatsM2(MARKERS, K, markereffectslist[[i]], P=Z,markermap, type=type, generation=generation)
      } else if (method==3){
        getstatsout<-getstatsM2(MARKERS, K, markereffectslist[[i]], P=Z)
      } else {
        getstatsout<-getstatsM1(MARKERS, K, markereffectslist[[i]], P=Z)
        
        }
      
      out<-c(out,list(-getstatsout[2], -getstatsout[3]))
    }
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
                    niterations,
                    type=type, generation=generation) 
  {
    
    Candidates<-1:nrow(allcombs)
    
    InitPop <- mclapply(1:npop, function(x) {
      return(sample(Candidates, ntoselect, replace=TRUE))
    }, mc.cores=mc.cores)
    
    InitPopFuncValues <- matrix(as.numeric(unlist(mclapply(InitPop,FUN = function(x) {
      GAobjfunc(x)}, mc.cores = mc.cores, mc.preschedule = T))), ncol= 2*length(markereffectslist)+1, byrow=TRUE)
    
    
    frontier3<- which(!emoa::is_dominated(t(InitPopFuncValues)))
    xy.f <- InitPopFuncValues[frontier3, ]
    colnames(xy.f)[1:length(colnamesforoutput)]<-colnamesforoutput
      
    if (plotiters){scatterplot3d::scatterplot3d(xy.f[,c(1,2,3)], highlight.3d=TRUE)}
    
    
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
                                                                }, mc.cores=mc.cores, mc.preschedule = T))), ncol= 2*length(markereffectslist)+1, byrow=TRUE)
      
      
      CurrentPop <- c(CurrentPop,ElitePop)
      notduplicated<-!duplicated(CurrentPop)
      CurrentPop<-CurrentPop[notduplicated]
      CurrentPopFuncValues<-rbind(CurrentPopFuncValues,ElitePopFuncValues)
      CurrentPopFuncValues<-CurrentPopFuncValues[notduplicated,]
      
      frontier3<- which(!emoa::is_dominated(t(CurrentPopFuncValues)))
      
      xy.f <- CurrentPopFuncValues[frontier3, ]
      
      
      #
      # Visualization.
      #
    colnames(xy.f)[1:length(colnamesforoutput)]<-colnamesforoutput
      
      if (plotiters){scatterplot3d::scatterplot3d(xy.f[,c(1,2,3)], highlight.3d=TRUE)}
      #points3d(CurrentPopFuncValues[,c(1,2,3)], col="blue")
      
      ElitePop <- mclapply(frontier3, FUN = function(x) {
        return(CurrentPop[[x]])
      }, mc.cores = mc.cores, mc.preschedule = T)
      ElitePopFuncValues <- CurrentPopFuncValues[frontier3,]
    }
    
    return(list(ElitePop, ElitePopFuncValues))
  }
  
  outGA<-GAfunc(ntoselect=nmates, npop=npopGA, mutprob=mutprob, 
                niterations=nitGA,
                type=type, generation=generation) 
  
  
  listofsols<-lapply(outGA[[1]], FUN=function(x){allcombs[x,]})
  nameslines<-rownames(MARKERS)
  
  if (!is.null(nameslines)){ if (length(nameslines)==nrow(MARKERS)){
    listofsols<-lapply(listofsols, FUN=function(x){
      return(cbind(nameslines[x[,1]],nameslines[x[,2]]))
    })
  }}
  outGA<-outGA[[2]]
  colnames(outGA)<-colnamesforoutput
  return(list(outGA, listofsols))
  
}



######################################################
## Group of functions that PopVar uses internally

par.position <- function(crossing.table, par.entries){ # Used when a crossing table is defined
  
  par.pos <- matrix(nrow = nrow(crossing.table), ncol = 2)
  crosses.possible <- matrix(nrow = nrow(crossing.table), ncol = 2)
  for(i in 1:nrow(crossing.table)){
    par1 <- as.character(crossing.table[i,1])
    par2 <- as.character(crossing.table[i,2])
    
    if(par1 %in% par.entries & par2 %in% par.entries){
      par.pos[i,1] <- which(par.entries == par1)
      par.pos[i,2] <- which(par.entries == par2)
      crosses.possible[i,1] <- par1
      crosses.possible[i,2] <- par2  
    }
  }
  
  
  for.dup <- paste(par.pos[,1], ".", par.pos[,2], sep=""); duplicated <- which(duplicated(for.dup))
  if(length(duplicated) > 0){
    par.pos <- par.pos[-duplicated, ]
    crosses.possible <- crosses.possible[-duplicated, ]
  }
  
  return(list(parent.positions=par.pos, crosses.possible=crosses.possible))
}

par.name <- function(crossing.mat, par.entries){ ## Used when all combinations of parents are crossed
  crosses.possible <- matrix(nrow = nrow(crossing.mat), ncol = 2)
  for(i in 1:nrow(crossing.mat)){
    crosses.possible[i,1] <- par.entries[crossing.mat[i,1]]
    crosses.possible[i,2] <- par.entries[crossing.mat[i,2]]
  }
  return(crosses.possible)
}

tails <- function(GEBVs, tail.p){ #Calculates means of tails; set tail.p to the proportion of the populaiton you want to take the mean of, default is 10%
  u.top <- mean(GEBVs[which(GEBVs >= quantile(GEBVs, 1-tail.p))], na.rm=TRUE)
  u.bot <- mean(GEBVs[which(GEBVs <= quantile(GEBVs, tail.p))], na.rm=TRUE)
  
  return(rbind(u.top, u.bot))
}




############################################


getstatsfromsim<-function (Markers, map, markereffectslist,crossing.table, tail.p = 0.2,  nSim = 25, simtype="riself")
{
  
  G.entries <- rownames(Markers)
  G.markers <- colnames(Markers)
  map.markers <- as.character(map[, 1])
  ntraits<-length(markereffectslist)
  
  crossing.mat <- par.position(crossing.table, par.entries = G.entries)$parent.position
  crosses.possible <- par.position(crossing.table, par.entries = G.entries)$crosses.possible
  
  
  
  ###########
  traits <- paste("T",1:length(markereffectslist),sep="_")
  name4out <- sample(10000:99999, 1)
  t.map <- t(map)
  rownames(t.map) <- NULL
  map4out <- cbind(c("pheno", "", ""), t.map)
  write.table(map4out, paste("map.tmp_", name4out, ".csv", 
                             sep = ""), row.names = F, col.names = F, sep = ",")
  options(warn = -1)
  read.map.out <- capture.output(read.map <- qtl::read.cross(format = "csv", 
                                                             crosstype = simtype, file = paste("map.tmp_", name4out, 
                                                                                                ".csv", sep = ""), na.strings = "NA"))
  
  unlink(paste("map.tmp_", name4out, ".csv", sep = ""))
  map_t1 <- qtl::pull.map(read.map)
  ########
  markereffectsmatrix<-as.matrix(Reduce("cbind", markereffectslist))
  
  par.BVs <- Markers%*%markereffectsmatrix 
  
  
  
  #################
  
  df.tmp <- data.frame(cbind(crosses.possible[, 1], crosses.possible[, 
                                                                     2], matrix(list(rep(NA, times = nSim)), nrow = nrow(crosses.possible), 
                                                                                ncol = (8 + 3 * (ntraits - 1)))))
  names(df.tmp)[1:10] <- c("Par1", "Par2", "midPar.Pheno", 
                           "midPar.GEBV", "pred.mu", "pred.mu_sd", "pred.varG", 
                           "pred.varG_sd", "mu.sp_low", "mu.sp_high")
  for (n in 1:ntraits){
    if (n == 1) 
      param.dfs <- list()
    param.dfs[[n]] <- as.matrix(df.tmp)
  }
  names(param.dfs) <- paste(traits, "_param.df", sep = "")
  
  ############
  p = 1
  M <- nrow(Markers)
  for (s in 1:nSim) {
    sim.pop <- qtl::sim.cross(map_t1, type = simtype, n.ind = M, 
                              model = NULL)
    qtl::write.cross(sim.pop, "csv", paste("sim.pop.tmp_", 
                                           name4out, sep = ""))
    pop.mat <- as.matrix(read.csv(paste("sim.pop.tmp_", name4out, 
                                        ".csv", sep = ""), header = T))[3:(M + 2), 2:(nrow(markereffectsmatrix) + 
                                                                                        1)]
    unlink(paste("sim.pop.tmp_", name4out, ".csv", sep = ""))
    for (z in 1:nrow(crossing.mat)) {
      pop.mat2 <- matrix(NA, nrow = nrow(pop.mat), ncol = ncol(pop.mat))
      par1 <- Markers[crossing.mat[z, 1], ]
      par2 <- Markers[crossing.mat[z, 2], ]
      for (r in 1:M) {
        pop.mat2[r, which(pop.mat[r, ] == "A")] <- par1[which(pop.mat[r, 
                                                                      ] == "A")]
        pop.mat2[r, which(pop.mat[r, ] == "B")] <- par2[which(pop.mat[r, 
                                                                      ] == "B")]
      }
      
      #mkr.has.0 <- apply(pop.mat2, 2, function(X) {
      #  return(length(which(X == 0)))
      #})
      #replace.0.mat <- rbind(which(mkr.has.0 != 0), mkr.has.0[which(mkr.has.0 != 
      #                                                                0)])
      #if (ncol(replace.0.mat) > 0) {
      #  for (b in 1:ncol(replace.0.mat)) {
      #    pop.mat2[which(pop.mat2[, replace.0.mat[1, 
      #                                            b]] == 0), replace.0.mat[1, b]] <- sample(c(1, 
      #                                                                                        -1), size = replace.0.mat[2, b], replace = T, 
      #                                                                                      prob = c(0.5, 0.5))
      #   }
      #  }
      prog_pred.mat <- pop.mat2 %*% markereffectsmatrix
      
      for (n in 1:ntraits) {
        if (s == 1) {
          if (ntraits > 1) 
            colnames(param.dfs[[n]])[11:(10 + 3 * (ntraits - 
                                                     1))] <- c(paste("low.resp_", traits[-n], 
                                                                     sep = ""), paste("high.resp_", traits[-n], 
                                                                                      sep = ""), paste("cor_w/_", traits[-n], 
                                                                                                       sep = ""))
          param.dfs[[n]][[z, "midPar.Pheno"]][s] <- 0
          param.dfs[[n]][[z, "midPar.GEBV"]][s] <- 0.5 * 
            (as.numeric(par.BVs[crossing.mat[z, 1], n]) + 
               as.numeric(par.BVs[crossing.mat[z, 2], 
                                  n]))
        }
        param.dfs[[n]][[z, "pred.mu"]][s] <- mean(prog_pred.mat[, 
                                                                n])
        param.dfs[[n]][[z, "pred.mu_sd"]][s] <- mean(prog_pred.mat[, 
                                                                   n])
        param.dfs[[n]][[z, "pred.varG"]][s] <- var(prog_pred.mat[, 
                                                                 n])
        param.dfs[[n]][[z, "pred.varG_sd"]][s] <- var(prog_pred.mat[, 
                                                                    n])
        param.dfs[[n]][[z, "mu.sp_low"]][s] <-  tails(prog_pred.mat[, 
                                                                    n], tail.p = tail.p)[2]
        param.dfs[[n]][[z, "mu.sp_high"]][s] <- tails(prog_pred.mat[, 
                                                                    n], tail.p = tail.p)[1]
        if (ntraits > 1) {
          index <- 1
          for (n2 in (1:ntraits)[-n]) {
            param.dfs[[n]][[z, 10 + index]][s] <- mean(prog_pred.mat[, 
                                                                     n2][which(prog_pred.mat[, n] <= quantile(prog_pred.mat[, 
                                                                                                                            n], probs = tail.p))], na.rm = T)
            param.dfs[[n]][[z, 10 + (ntraits - 1) + index]][s] <- mean(prog_pred.mat[, 
                                                                                     n2][which(prog_pred.mat[, n] >= quantile(prog_pred.mat[, 
                                                                                                                                            n], probs = 1 - tail.p))], na.rm = T)
            param.dfs[[n]][[z, 10 + 2 * (ntraits - 1) + 
                              index]][s] <- cor(prog_pred.mat[, n], prog_pred.mat[, 
                                                                                  n2], use = "complete.obs")
            index <- index + 1
          }
        }
      }
      p <- p + 1
    }
  }
  
  
  
  
  
  ######
  preds.per.sim <- param.dfs
  for (n in 1:ntraits) {
    col.names <- colnames(param.dfs[[n]])
    for (c in 3:length(col.names)) {
      name.tmp <- col.names[c]
      if (name.tmp %in% c("pred.mu_sd", "pred.varG_sd")) 
        param.dfs[[n]][, c] <- sapply(param.dfs[[n]][, 
                                                     c], FUN = sd, na.rm = T)
      if (!name.tmp %in% c("pred.mu_sd", "pred.varG_sd")) 
        param.dfs[[n]][, c] <- sapply(param.dfs[[n]][, 
                                                     c], FUN = mean, na.rm = T)
    }
  }
  
  return(list(predictions = param.dfs, pop.mat2=pop.mat2))
  
}



################################################


getGaSolutionsFrontierMultiTraitSimcross<-function(Markers, Markers2=NULL,K, map,markereffectslist,nmates=NULL,  nSim = 5,npopGA, nitGA, mutprob, mc.cores, noself=F,
                                                   simtype="riself", plotiters=F){
  ##############################list of parents optimization
  tail.p = .00001 
  
  if(is.null(names(markereffectslist))){
    
    names(markereffectslist)<-paste("T_", 1:length(markereffectslist),sep="")
  }
  
  
  colnamesforoutput<-c("I",paste(rep(c("-G_", "-U_"),length(names(markereffectslist))),rep(names(markereffectslist), each=length(markereffectslist)),sep=""))  
  
  if (!is.null(Markers2)){
    N1<-nrow(Markers)
    N2<-nrow(Markers2)
    
    MARKERS<-rbind(Markers,Markers2)
  } else {
    MARKERS=Markers
  }
  
  N=nrow(MARKERS)
  
  if (is.null(Markers2)){
    if (is.null(nmates)){nmates=N}
    allcombs<-matrix(nrow=0,ncol=2)
    ij=1
    if (noself){ 
      for (i in 1:(N-1)){for(j in (i+1):N){allcombs<-rbind(allcombs,c(i,j));ij+ij+1}}
    } else {
      for (i in 1:N){for(j in i:N){allcombs<-rbind(allcombs,c(i,j));ij+ij+1}}
    }
    
  } else {
    if (is.null(nmates)){nmates=N1}
    allcombs<-matrix(nrow=0,ncol=2)
    ij=1
    for (i in 1:N1){for(j in (N1+1):(N1+N2)){allcombs<-rbind(allcombs,c(i,j));ij+ij+1}}
  }
  
  GAobjfunc<-function(indexesinallcombs){
    Z<-matrix(0,nrow=nmates, ncol=N)
    for (i in 1:length(indexesinallcombs)){
      Z[i,allcombs[indexesinallcombs[i],]]<-1/2
      if (allcombs[indexesinallcombs[i],][1]==allcombs[indexesinallcombs[i],][2]){Z[i,allcombs[indexesinallcombs[i],]]<-1}
    }
    
    suppressWarnings(suppressMessages(out0<-getstatsfromsim(Markers=MARKERS, map=map, markereffectslist=markereffectslist,crossing.table=allcombs[indexesinallcombs,], tail.p = tail.p ,  nSim = nSim, simtype=simtype)))
    
    getstatsout<-getstatsM1(MARKERS, K, markereffectslist[[1]], P=Z)
    
    out1<-as.matrix(Reduce("cbind",out0$predictions))
    out1<-out1[,c(1,2,which(!colnames(out1)%in%c("Par1", "Par2", "midPar.Pheno")))]
    out1<-out1[match(paste(allcombs[indexesinallcombs,1],allcombs[indexesinallcombs,2],sep=""),paste(out1[,1], out1[,2],sep="")),]
    
    out1<-out1[,-c(1,2)]
    namesout1<-colnames(out1)
    
    namesout1in<-which(namesout1%in%c("midPar.GEBV", "pred.mu","pred.varG")) 
    out1<-matrix(as.numeric(out1), nrow=nrow(out1))
    out2<-apply(out1[,namesout1in],2,function(x){return(sum(x))})
    names(out2)<-paste(namesout1[namesout1in],rep(1:length(markereffectslist), each=3),sep="_")

    out3<-getstatsout[1]
    
    #############error here 
    for (i in 1:length(names(out2))){
      out3<-c(out3,-out2[i])
    }
    out3<-out3[c(T,((1:(length(out3)-1))%%3) %in% c(2,0))]
    
    return(out3)
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
      return(sample(Candidates, ntoselect, replace=TRUE))
    }, mc.cores=mc.cores)
    
    InitPopFuncValues <- matrix(as.numeric(unlist(mclapply(InitPop,FUN = function(x) {
      GAobjfunc(x)}, mc.cores = mc.cores, mc.preschedule = T))), ncol= 1+2*length(markereffectslist), byrow=TRUE)
    
    
    frontier3<- which(!emoa::is_dominated(t(InitPopFuncValues)))
    xy.f <- InitPopFuncValues[frontier3, ]
    
   
    colnames(xy.f)<-colnamesforoutput
    
    if (plotiters){scatterplot3d::scatterplot3d(xy.f[,c(1,2,3)], highlight.3d=TRUE)}
    
    
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
                                                                }, mc.cores=mc.cores, mc.preschedule = T))), ncol= 1+2*length(markereffectslist), byrow=TRUE)
      
      
      CurrentPop <- c(CurrentPop,ElitePop)
      notduplicated<-!duplicated(CurrentPop)
      CurrentPop<-CurrentPop[notduplicated]
      CurrentPopFuncValues<-rbind(CurrentPopFuncValues,ElitePopFuncValues)
      CurrentPopFuncValues<-CurrentPopFuncValues[notduplicated,]
      
      frontier3<- which(!emoa::is_dominated(t(CurrentPopFuncValues)))
      
      xy.f <- CurrentPopFuncValues[frontier3, ]
    colnames(xy.f)<-colnamesforoutput
      
      
      #
      # Visualization.
      #
      if (plotiters){scatterplot3d::scatterplot3d(xy.f[,c(1,2,3)], highlight.3d=TRUE)}
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
  nameslines<-rownames(MARKERS)
  
  if (!is.null(nameslines)){ if (length(nameslines)==nrow(MARKERS)){
    listofsols<-lapply(listofsols, FUN=function(x){
      return(cbind(nameslines[x[,1]],nameslines[x[,2]]))
    })
  }}
  outGA<-outGA[[2]]
  colnames(outGA)<-colnamesforoutput
  
  return(list(outGA, listofsols))
  
}


############################################

pairs3d<-function(X, cols=c(1,2,3,4),mfrow=c(2,2)){
  colnamesX<-colnames(X)
  if(is.null(colnamesX)){
    colnames(X)<-colnamesX<-1:ncol(X)
  }
  allcombs<-combn(colnamesX, 3)
  allcombsnum<-combn(1:length(colnamesX), 3)
  
  lengthplot<- ncol(allcombsnum)
  nrowplot<-lengthplot/length(colnamesX)
  ncolplot<-lengthplot/nrowplot
  par(mfrow=mfrow)
  for (i in 1:ncol(allcombsnum)){
    scatterplot3d::scatterplot3d(X[,allcombsnum[,i]])
  }
  
  par(mfrow=c(1,1))
  
  return()
}

############################################

plotGM<-function(GMsols, type="3D",plotly=FALSE, idealsol=NULL, traitnum=1){
  
 
  
  if (type=="3D"){
    if (is.null(idealsol)){
      gasolsstandardized<-scale(GMsols[[1]])
      Dsols<-as.matrix(dist(rbind(c(min(gasolsstandardized[,1]),min(gasolsstandardized[,2*(traitnum-1)+2]), min(gasolsstandardized[,2*(traitnum-1)+3])),gasolsstandardized[,c(1,2*(traitnum-1)+2,2*(traitnum-1)+3)])))
      disttobest<-Dsols[-1,1] 
      colvar<-disttobest^2 
      colvar=(colvar-min(colvar))/(max(colvar)-min(colvar))
      textforplot = colnames(GMsols[[1]])
      textforplot<-apply(GMsols[[1]],1,function(x){paste(c(paste(textforplot,x)),sep="",collapse="\n")})
      textforplot<-paste(textforplot,paste("\n obs", 1:nrow(GMsols[[1]]),sep=""))

      
      # Rank variable for colour assignment
      if (plotly){
      dforder = findInterval(exp(-colvar), sort(exp(-colvar)))
      p1 <- plotly::plot_ly(x = GMsols[[1]][,1],y=GMsols[[1]][,2*(traitnum-1)+2], z=GMsols[[1]][,2*(traitnum-1)+3],mode="markers",
                    marker = list(color = exp(-colvar), colorscale = c('#FFE1A1', '#EADDDA'), showscale = FALSE),
                    text = textforplot,
                    hoverinfo="text") %>%
        plotly::add_markers() %>%
        plotly::layout(scene = list(xaxis = list(title = 'I'),
                            yaxis = list(title = '-G'),
                            zaxis = list(title = '-U'))
        ) %>%
        plotly::layout(showlegend = FALSE)
      
      
      p1%>%plotly::offline(height=400)
      } else {tempdata<-data.frame(GMsols[[1]][,1],GMsols[[1]][,2*(traitnum-1)+2], GMsols[[1]][,2*(traitnum-1)+3])
      scatterplot3d::scatterplot3d(tempdata, highlight.3d=T,
                                   xlab="I", ylab="-G", zlab="-U")
              }
    } else{
      gasolsstandardized<-GMsols[[1]]
      Dsols<-as.matrix(dist(rbind(idealsol,gasolsstandardized)))
      disttobest<-Dsols[-1,1]
      colvar<-disttobest^2
      colvar=(colvar-min(colvar))/(max(colvar)-min(colvar))
      
      
      textforplot = colnames(GMsols[[1]])
      textforplot<-apply(GMsols[[1]],1,function(x){paste(c(paste(textforplot,x)),sep="",collapse="\n")})
      textforplot<-paste(textforplot,paste("\n obs", 1:nrow(GMsols[[1]]),sep=""))
      if (plotly){
      # Rank variable for colour assignment
      dforder = findInterval(exp(-colvar), sort(exp(-colvar)))
      p1 <- plotly::plot_ly(x = GMsols[[1]][,1],y=GMsols[[1]][,2*(traitnum-1)+2], z=GMsols[[1]][,2*(traitnum-1)+3],mode="markers",
                    marker = list(color = exp(-colvar), colorscale = c('#FFE1A1', '#EADDDA'), showscale = FALSE),
                    text = textforplot,
                    hoverinfo="text") %>%
        plotly::add_markers() %>%
        plotly::layout(scene = list(xaxis = list(title = 'I'),
                            yaxis = list(title = '-G'),
                            zaxis = list(title = '-U'))
        ) %>%
        plotly::layout(showlegend = FALSE)
      
      
      p1%>%plotly::offline(height=400)
      
      } else {
        tempdata<-data.frame(GMsols[[1]][,1],GMsols[[1]][,2*(traitnum-1)+2], GMsols[[1]][,2*(traitnum-1)+3])
        scatterplot3d::scatterplot3d(tempdata,  highlight.3d=T,
                                     xlab="I", ylab="-G", zlab="-U")
        }
    }
  }
  if (type=="SOM"){
    gasolsstandardized<-scale(GMsols[[1]])
    gasolsstandardized2<-apply(gasolsstandardized,2,function(x){(x-min(x))/(max(x)-min(x))})
    gasolsstandardized2[,2*(traitnum-1)+2]<--gasolsstandardized2[,2*(traitnum-1)+2]
    gasolsstandardized2[,2*(traitnum-1)+3]<--gasolsstandardized2[,2*(traitnum-1)+3]
    colnames(gasolsstandardized2)<-unlist(lapply(colnames(GMsols[[1]]),function(x){gsub("-","",x)}))
    
    GMsols.SOM1 <- kohonen::som(gasolsstandardized2[,c(1,2*(traitnum-1)+2,2*(traitnum-1)+3)], grid = somgrid(10, 10, "rectangular"))
    
    par(mfrow = c(1, 2))
    plot(GMsols.SOM1, type = "mapping", pchs = 0, main = "Mapping Type SOM",labels = 1:nrow(gasolsstandardized2),cex=.5 )
    plot(GMsols.SOM1, main = "SOM values",palette.name = topo.colors)
    par(mfrow = c(1, 1))
  }
  if (type=="SOM2"){
    gasolsstandardized<-scale(GMsols[[1]])
    gasolsstandardized2<-apply(gasolsstandardized,2,function(x){(x-min(x))/(max(x)-min(x))})
    gasolsstandardized2[,2*(traitnum-1)+2]<--gasolsstandardized2[,2*(traitnum-1)+2]
    gasolsstandardized2[,2*(traitnum-1)+3]<--gasolsstandardized2[,2*(traitnum-1)+3]
    colnames(gasolsstandardized2)<-unlist(lapply(colnames(GMsols[[1]]),function(x){gsub("-","",x)}))
    #GMsols.SOM1 <- som(gasolsstandardized2, grid = somgrid(6, 4, "rectangular"))
    #
    default.paramSOM <- initSOM(dimension=c(5,5))
    tempdat<-gasolsstandardized2[,c(1,2*(traitnum-1)+2,2*(traitnum-1)+3)]
    colnames(tempdat)<-c("I","G","U")
    som.lemis <-SOMbrero::trainSOM(tempdat , maxit=500, radius.type="letremy")
    par(mfrow=c(2,2))
    plot(som.lemis, what="prototypes", type="3d", variable=1, print.title=TRUE, 
         main="")
    plot(som.lemis, what="prototypes", type="3d", variable=2, print.title=TRUE, 
         main="")
    plot(som.lemis, what="prototypes", type="3d", variable=3, print.title=TRUE, 
         main="")
    plot(som.lemis,  what="obs", type="radar", print.title=F, 
         main="")
    par(mfrow = c(1, 1))
  }
}



