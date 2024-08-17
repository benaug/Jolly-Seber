e2dist  <-  function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.JS.SCR.SexPopDy <- function(lambda.y1.M=NA,lambda.y1.F=NA,n.year=NA,
                                gamma.sex=NA,phi.sex=NA,
                   p0.sex=NA,sigma.sex=NA,X=NA,buff=buff,K=NA,
                   p.obs.sex=NA,sigma.move.sex=NULL){
  #Population dynamics
  N <- N.M <- N.F <- rep(NA,n.year)
  N.recruit.M <- N.survive.M <- ER.M <- rep(NA,n.year-1)
  N.recruit.F <- N.survive.F <- ER.F <- rep(NA,n.year-1)
  N.M[1] <- rpois(1,lambda.y1.M)
  N.F[1] <- rpois(1,lambda.y1.F)
  N[1] <- N.M[1]+N.F[1]
  
  #Easiest to increase dimension of z as we simulate bc size not known in advance.
  z <- matrix(0,N[1],n.year)
  z[1:N[1],1] <- 1
  sex <- c(rep(1,N.M[1]),rep(2,N.F[1])) #derived variable
  phi <- matrix(NA,N[1],n.year-1)
  for(g in 2:n.year){
    #Simulate recruits
    ER.M[g-1] <- N[g-1]*gamma.sex[1]
    ER.F[g-1] <- N[g-1]*gamma.sex[2]
    N.recruit.M[g-1] <- rpois(1,ER.M[g-1])
    N.recruit.F[g-1] <- rpois(1,ER.F[g-1])
    #add recruits to z
    z.dim.old <- nrow(z)
    z <- rbind(z,matrix(0,nrow=N.recruit.M[g-1]+N.recruit.F[g-1],ncol=n.year))
    z[(z.dim.old+1):(z.dim.old+N.recruit.M[g-1]+N.recruit.F[g-1]),g] <- 1
    #record sexes
    sex <- c(sex,rep(1,N.recruit.M[g-1]),rep(2,N.recruit.F[g-1]))
    
    #Simulate survival
    phi <- rbind(phi,matrix(NA,nrow=N.recruit.M[g-1]+N.recruit.F[g-1],ncol=n.year-1))
    phi[,g-1] <- phi.sex[sex]
    
    idx <- which(z[,g-1]==1)
    z[idx,g] <- rbinom(length(idx),1,phi[idx,g-1])
    N.survive.M[g-1] <- sum(z[,g-1]==1&z[,g]==1&sex==1)
    N.survive.F[g-1] <- sum(z[,g-1]==1&z[,g]==1&sex==2)
    N.M[g] <- N.recruit.M[g-1]+N.survive.M[g-1]
    N.F[g] <- N.recruit.F[g-1]+N.survive.F[g-1]
    N[g] <- N.M[g]+N.F[g]
  }
  N.recruit <- N.recruit.M+N.recruit.F
  N.survive <- N.survive.M+N.survive.F
  
  if(any(N.recruit.M+N.survive.M!=N.M[2:n.year]))stop("Simulation bug")
  if(any(N.recruit.F+N.survive.F!=N.F[2:n.year]))stop("Simulation bug")
  if(any(colSums(z)!=N))stop("Simulation bug")
  
  #detection
  #get maximal x and y extent across yearly grids plus buffer
  xlim <- c(max(unlist(lapply(X,function(x){min(x[,1])}))),max(unlist(lapply(X,function(x){max(x[,1])})))) + c(-buff,buff)
  ylim <- c(max(unlist(lapply(X,function(x){min(x[,2])}))),max(unlist(lapply(X,function(x){max(x[,2])})))) + c(-buff,buff)
  J <- unlist(lapply(X,nrow)) #extract number of traps per year
  J.max <- max(J)
  
  #simulate activity centers - fixed through time
  N.super <- nrow(z)
  library(truncnorm)
  if(!is.null(sigma.move.sex)){
    print("simulating mobile ACs")
    s <- array(NA,dim=c(N.super,n.year,2))
    s[,1,] <- cbind(runif(N.super, xlim[1],xlim[2]), runif(N.super,ylim[1],ylim[2]))
    for(g in 2:n.year){
      for(i in 1:N.super){
        s[i,g,1] <- rtruncnorm(1,s[i,g-1,1],sd=sigma.move.sex[sex[i]],a=xlim[1],b=xlim[2])
        s[i,g,2] <- rtruncnorm(1,s[i,g-1,2],sd=sigma.move.sex[sex[i]],a=ylim[1],b=ylim[2])
      }
    }
  }else{
    print("simulating fixed ACs (provide sigma.move.sex for mobile)")
    s <- cbind(runif(N.super, xlim[1],xlim[2]), runif(N.super,ylim[1],ylim[2]))
  }
  
  pd <- y <- array(0,dim=c(N.super,n.year,J.max))
  
  for(g in 1:n.year){
    if(!is.null(sigma.move.sex)){
      D <- e2dist(s[,g,],X[[g]])
    }else{
      D <- e2dist(s,X[[g]])
    }
    pd[,g,1:J[g]] <- p0.sex[sex]*exp(-D*D/(2*sigma.sex[sex]^2))
    for(i in 1:N.super){
      if(z[i,g]==1){
          y[i,g,1:J[g]] <- rbinom(J[g],K[g],pd[i,g,1:J[g]])
      }
    }
  }
  
  #store true data for model buildling/debugging
  truth <- list(y=y,cov=cov,N=N,N.recruit=N.recruit,N.survive=N.survive,
             N.M=N.M,N.F=N.F,N.recruit.M=N.recruit.M,N.recruit.F=N.recruit.F,
             N.survive.M=N.survive.M,N.survive.F=N.survive.F,
             sex=sex,z=z,s=s)
  
  #discard undetected individuals
  keep.idx <- which(rowSums(y)>0)
  y <- y[keep.idx,,]
  sex <- sex[keep.idx]
  
  #discard unobserved sexes
  is.missing <- rbinom(length(sex),1,1-p.obs.sex)
  sex[which(is.missing==1)] <- NA
  
  return(list(y=y,sex=sex,
              N=N,N.recruit=N.recruit,N.survive=N.survive,
              N.M=N.M,N.recruit.M=N.recruit.M,N.survive.M=N.survive.M,
              N.F=N.F,N.recruit.F=N.recruit.F,N.survive.F=N.survive.F,
              X=X,K=K,
              xlim=xlim,ylim=ylim,truth=truth))
}
