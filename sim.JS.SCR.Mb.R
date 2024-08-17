e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.JS.SCR.Mb <- function(lambda.y1=NA,gamma=NA,n.year=NA,
                       beta0.phi=NA,beta1.phi=NA,
                   p01=NA,p02=NA,sigma=NA,X=NA,buff=buff,K=NA){
  #Population dynamics
  N <- rep(NA,n.year)
  N.recruit <- N.survive <- ER <- rep(NA,n.year-1)
  N[1] <- rpois(1,lambda.y1)
  
  #Easiest to increase dimension of z as we simulate bc size not known in advance.
  z <- matrix(0,N[1],n.year)
  z[1:N[1],1] <- 1
  cov <- rnorm(N[1],0,1) #simulate ind survival covariate for 1st year guys
  phi <- matrix(NA,N[1],n.year-1)
  for(g in 2:n.year){
    #Simulate recruits
    ER[g-1] <- N[g-1]*gamma[g-1]
    N.recruit[g-1] <- rpois(1,ER[g-1])
    #add recruits to z
    z.dim.old <- length(cov)
    z <- rbind(z,matrix(0,nrow=N.recruit[g-1],ncol=n.year))
    z[(z.dim.old+1):(z.dim.old+N.recruit[g-1]),g] <- 1
    cov <- c(cov,rep(NA,N.recruit[g-1]))
    cov[(z.dim.old+1):(z.dim.old+N.recruit[g-1])] <- rnorm(N.recruit[g-1],0,1) #simulate survival cov values for new recruits
    
    #Simulate survival
    phi <- rbind(phi,matrix(NA,nrow=N.recruit[g-1],ncol=n.year-1))
    phi[,g-1] <- plogis(beta0.phi+cov*beta1.phi)
    
    idx <- which(z[,g-1]==1)
    z[idx,g] <- rbinom(length(idx),1,phi[idx,g-1])
    N.survive[g-1] <- sum(z[,g-1]==1&z[,g]==1)
    N[g] <- N.recruit[g-1]+N.survive[g-1]
  }
  
  if(any(N.recruit+N.survive!=N[2:n.year]))stop("Simulation bug")
  if(any(colSums(z)!=N))stop("Simulation bug")
  
  #plot to see if sim values realistic
  hist(phi,main="Distribution of Individual Phi")
  
  #detection
  #get maximal x and y extent across yearly grids plus buffer
  xlim <- c(max(unlist(lapply(X,function(x){min(x[,1])}))),max(unlist(lapply(X,function(x){max(x[,1])})))) + c(-buff,buff)
  ylim <- c(max(unlist(lapply(X,function(x){min(x[,2])}))),max(unlist(lapply(X,function(x){max(x[,2])})))) + c(-buff,buff)
  J <- unlist(lapply(X,nrow)) #extract number of traps per year
  J.max <- max(J)
  
  #simulate activity centers - fixed through time
  N.super <- nrow(z)
  s<- cbind(runif(N.super, xlim[1],xlim[2]), runif(N.super,ylim[1],ylim[2]))
  K.max <- max(K)
  y <- array(0,dim=c(N.super,n.year,J.max,K.max))
  pd1 <- pd2 <- array(0,dim=c(N.super,n.year,J.max))
  
  for(g in 1:n.year){
    D<- e2dist(s,X[[g]])
    pd1[,g,1:J[g]]<- p01[g]*exp(-D*D/(2*sigma[g]*sigma[g]))
    pd2[,g,1:J[g]]<- p02[g]*exp(-D*D/(2*sigma[g]*sigma[g]))
    for(i in 1:N.super){
      if(z[i,g]==1){
        state <- matrix(0,J[g],K[g]) #one individual's trap by occasion capture state history
        for(j in 1:J[g]){
          for(k in 1:K[g]){
            if(state[j,k]==0){
              y[i,g,j,k] <- rbinom(1,1,pd1[i,g,j]) #pd doesn't vary by occasion
            }else{
              y[i,g,j,k] <- rbinom(1,1,pd2[i,g,j])
            }
            if((y[i,g,j,k]==1)&(k<K[g])){
              state[j,(k+1):K[g]] <- 1
            }
          }
        }
      }
    }
  }
  
  #store true data for model buildling/debugging
  truth <- list(y=y,cov=cov,N=N,N.recruit=N.recruit,N.survive=N.survive,z=z,s=s)
  
  #discard undetected individuals
  keep.idx <- which(rowSums(y)>0)
  y <- y[keep.idx,,,]
  cov <- cov[keep.idx]
  return(list(y=y,cov=cov,N=N,N.recruit=N.recruit,N.survive=N.survive,X=X,K=K,
              xlim=xlim,ylim=ylim,truth=truth))
}
