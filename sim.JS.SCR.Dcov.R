e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.JS.SCR.Dcov <- function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,
                            gamma=NA,n.year=NA,beta0.phi=NA,beta1.phi=NA,
                   p0=NA,sigma=NA,X=NA,buff=buff,K=NA,xlim=NA,ylim=NA,res=NA){
  #Population dynamics
  N <- rep(NA,n.year)
  N.recruit <- N.survive <- ER <- rep(NA,n.year-1)
  #get expected N in year 1 from D.cov parameters
  cellArea <- res^2
  lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
  lambda.y1 <- sum(lambda.cell)
  N[1] <- rpois(1,lambda.y1)

  #recreate some Dcov things so we can pass fewer arguments into this function
  x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
  y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
  dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
  cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
  n.cells <- nrow(dSS)
  n.cells.x <- length(x.vals)
  n.cells.y <- length(y.vals)

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

  #simulate activity centers - fixed through time
  N.super <- nrow(z)
  # simulate a population of activity centers
  pi.cell <- lambda.cell/sum(lambda.cell)
  s.cell <- sample(1:n.cells,N.super,prob=pi.cell,replace=TRUE)
  #distribute activity centers uniformly inside cells
  s <- matrix(NA,nrow=N.super,ncol=2)
  for(i in 1:N.super){
    s.xlim <- dSS[s.cell[i],1] + c(-res,res)/2
    s.ylim <- dSS[s.cell[i],2] + c(-res,res)/2
    s[i,1] <- runif(1,s.xlim[1],s.xlim[2])
    s[i,2] <- runif(1,s.ylim[1],s.ylim[2])
  }

  pd <- y <- array(0,dim=c(N.super,n.year,J.max))
  for(g in 1:n.year){
    D <- e2dist(s,X[[g]])
    pd[,g,1:J[g]] <- p0[g]*exp(-D*D/(2*sigma[g]*sigma[g]))
    for(i in 1:N.super){
      if(z[i,g]==1){
          y[i,g,1:J[g]] <- rbinom(J[g],K[g],pd[i,g,1:J[g]])
      }
    }
  }

  #store true data for model building/debugging
  truth <- list(y=y,cov=cov,N=N,N.recruit=N.recruit,N.survive=N.survive,z=z,s=s)

  #discard undetected individuals
  keep.idx <- which(rowSums(y)>0)
  y <- y[keep.idx,,]
  cov <- cov[keep.idx]
  return(list(y=y,cov=cov,N=N,N.recruit=N.recruit,N.survive=N.survive,X=X,K=K,
              xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
              n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,
              D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N,lambda.y1=lambda.y1,
              truth=truth))
}
