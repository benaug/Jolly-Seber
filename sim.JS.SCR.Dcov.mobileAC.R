e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.JS.SCR.Dcov.mobileAC <- function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,
                            gamma=NA,n.year=NA,beta0.phi=NA,beta1.phi=NA,
                   p0=NA,sigma=NA,sigma.move=sigma.move,rsf.beta=rsf.beta,
                   X=NA,buff=buff,K=NA,xlim=NA,ylim=NA,res=NA){
  #Population dynamics
  N <- rep(NA,n.year)
  N.recruit <- N.survive <- ER <- rep(NA,n.year-1)
  #get expected N in year 1 from D.cov parameters
  cellArea <- res^2
  lambda.cell <- exp(D.beta0 + D.beta1*D.cov)*cellArea
  lambda.y1 <- sum(lambda.cell)
  N[1] <- rpois(1,lambda.y1)

  #recreate some Dcov things so we can pass fewer arguments into this function
  x.vals <- seq(xlim[1],xlim[2],by=res)
  y.vals <- seq(ylim[1],ylim[2],by=res)
  dSS <- as.matrix(expand.grid(x.vals,y.vals)) + res/2 #add res/2 to get cell centroids
  #remove extra cells created outside xlim and ylim
  rem.idx <- which(dSS[,1]>xlim[2]|dSS[,2]>ylim[2])
  dSS <- dSS[-rem.idx,]
  cells <- matrix(1:nrow(dSS),nrow=length(x.vals)-1,ncol=length(y.vals)-1)
  n.cells <- nrow(dSS)
  n.cells.x <- length(x.vals) - 1
  n.cells.y <- length(y.vals) - 1

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
  J <- unlist(lapply(X,nrow)) #extract number of traps per year
  J.max <- max(J)
  N.super <- nrow(z)
  
  # simulate a population of activity centers for year 1 proportional to D.cov
  library(truncnorm)
  pi.cell <- lambda.cell/sum(lambda.cell)
  #zero out non-habitat
  pi.cell[InSS==0] <- 0
  s.cell <- matrix(NA,N.super,n.year)
  s.cell[,1] <- sample(1:n.cells,N.super,prob=pi.cell,replace=TRUE)
  #distribute activity centers uniformly inside cells
  s <- array(NA,dim=c(N.super,n.year,2))
  for(i in 1:N.super){
    tmp <- which(cells==s.cell[i,1],arr.ind=TRUE) #x and y number
    s[i,1,1] <- runif(1,x.vals[tmp[1]],x.vals[tmp[1]+1])
    s[i,1,2] <- runif(1,y.vals[tmp[2]],y.vals[tmp[2]+1])
  }
  #subsequent years
  avail.dist <- use.dist <- array(NA,dim=c(N.super,n.year-1,n.cells))
  rsf <- exp(rsf.beta*D.cov)
  for(g in 2:n.year){
    for(i in 1:N.super){
      avail.dist[i,g-1,] <- getAvail(s=s[i,g-1,1:2],sigma=sigma.move,res=res,x.vals=x.vals,
                                 y.vals=y.vals,n.cells.x=n.cells.x,n.cells.y=n.cells.y)
      use.dist[i,g-1,] <- rsf*avail.dist[i,g-1,]
      use.dist[i,g-1,] <- use.dist[i,g-1,]/sum(use.dist[i,g-1,])
      #move AC - select new cell
      s.cell[i,g] <- sample(1:n.cells,1,replace=TRUE,prob=use.dist[i,g-1,])
      #choose location inside cell
      s.xlim <- dSS[s.cell[i,g],1] + c(-res,res)/2
      s.ylim <- dSS[s.cell[i,g],2] + c(-res,res)/2
      #choose new location inside cell
      s[i,g,1] <- rtruncnorm(1,a=s.xlim[1],b=s.xlim[2],mean=s[i,g-1,1],sd=sigma.move)
      s[i,g,2] <- rtruncnorm(1,a=s.ylim[1],b=s.ylim[2],mean=s[i,g-1,2],sd=sigma.move)
    }
  }
  pd <- y <- array(0,dim=c(N.super,n.year,J.max))
  for(g in 1:n.year){
    D <- e2dist(s[,g,],X[[g]])
    pd[,g,1:J[g]]<- p0[g]*exp(-D*D/(2*sigma[g]*sigma[g]))
    for(i in 1:N.super){
      if(z[i,g]==1){
        y[i,g,1:J[g]] <- rbinom(J[g],K[g],pd[i,g,1:J[g]])
      }
    }
  }
  
  #expected proportion of realized N in cell 
  pi.cell <- array(NA,dim=c(n.year,n.cells))
  pi.cell[1,] <- lambda.cell/sum(lambda.cell)
  for(g in 2:n.year){
    pi.cell[g,] <- colSums(use.dist[z[,g]==1,g-1,])
  }

  #store true data for model building/debugging
  truth <- list(y=y,cov=cov,N=N,N.recruit=N.recruit,N.survive=N.survive,z=z,s=s,pi.cell=pi.cell)

  #discard undetected individuals
  keep.idx <- which(rowSums(y)>0)
  y <- y[keep.idx,,]
  cov <- cov[keep.idx]
  return(list(y=y,cov=cov,N=N,N.recruit=N.recruit,N.survive=N.survive,X=X,K=K,
              xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
              n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,s=s,
              D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N,lambda.y1=lambda.y1,
              truth=truth))
}
