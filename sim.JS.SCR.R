e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.JS.SCR <- function(lambda.y1=NA,gamma=NA,n.year=NA,
                       beta0.phi=NA,beta1.phi=NA,
                   p0=NA,sigma=NA,X=NA,buff=buff,K=NA){
  #Population dynamics
  N=N.recruit=N.survive=ER=rep(NA,n.year)
  N[1]=rpois(1,lambda.y1)
  
  #Easiest to increase dimension of z as we simulate bc size not known in advance.
  z=matrix(0,N[1],n.year)
  z[1:N[1],1]=1
  cov=rnorm(N[1],0,1) #simulate ind survival covariate for 1st year guys
  phi=matrix(NA,N[1],n.year-1) #to store individual survival covariate
  for(g in 2:n.year){
    #Simulate recruits
    ER[g]=N[g-1]*gamma[g]
    N.recruit[g]=rpois(1,ER[g])
    #add recruits to z
    z.dim.old=length(cov)
    z=rbind(z,matrix(0,nrow=N.recruit[g],ncol=n.year))
    z[(z.dim.old+1):(z.dim.old+N.recruit[g]),g]=1
    cov=c(cov,rep(NA,N.recruit[g]))
    cov[(z.dim.old+1):(z.dim.old+N.recruit[g])]=rnorm(N.recruit[g],0,1) #simulate survival cov values for new recruits
    
    #Simulate survival
    phi=rbind(phi,matrix(NA,nrow=N.recruit[g],ncol=n.year-1))
    phi[,g-1]=plogis(beta0.phi+cov*beta1.phi) #no year effects here, just replicating across years
    
    idx=which(z[,g-1]==1)
    z[idx,g]=rbinom(length(idx),1,phi[idx,g-1])
    N.survive[g]=sum(z[,g-1]==1&z[,g]==1)
    N[g]=N.recruit[g]+N.survive[g]
  }
  
  if(any(N.recruit[2:n.year]+N.survive[2:n.year]!=N[2:n.year]))stop("Simulation bug")
  if(any(colSums(z)!=N))stop("Simulation bug")
  
  #plot to see if sim values realistic
  hist(phi,main="dist of phi across individuals")
  
  birth.year=apply(z,1,function(x){match(1,x)})
  death.year=10-apply(z,1,function(x){match(1,rev(x))})+1
  lifetimes=death.year-birth.year
  
  #detection
  #get maximal x and y extent across yearly grids plus buffer
  xlim = c(max(unlist(lapply(X,function(x){min(x[,1])}))),max(unlist(lapply(X,function(x){max(x[,1])})))) + c(-buff,buff)
  ylim = c(max(unlist(lapply(X,function(x){min(x[,2])}))),max(unlist(lapply(X,function(x){max(x[,2])})))) + c(-buff,buff)
  J = unlist(lapply(X,nrow)) #extract number of traps per year
  J.max=max(J)
  
  #simulate activity centers - fixed through time
  N.super=nrow(z)
  s<- cbind(runif(N.super, xlim[1],xlim[2]), runif(N.super,ylim[1],ylim[2]))
  pd <- y <- array(0,dim=c(N.super,n.year,J.max))
  
  for(g in 1:n.year){
    D<- e2dist(s,X[[g]])
    pd[,g,1:J[g]]<- p0[g]*exp(-D*D/(2*sigma[g]*sigma[g]))
    for(i in 1:N.super){
      if(z[i,g]==1){
          y[i,g,1:J[g]] <- rbinom(J[g],K[g],pd[i,g,1:J[g]])
      }
    }
  }
  
  #store true data for model buildling/debugging
  truth=list(y=y,cov=cov,N=N,N.recruit=N.recruit,N.survive=N.survive,z=z,s=s)
  
  #discard undetected individuals
  keep.idx=which(rowSums(y)>0)
  y=y[keep.idx,,]
  cov=cov[keep.idx]
  return(list(y=y,cov=cov,N=N,N.recruit=N.recruit,N.survive=N.survive,X=X,K=K,
              xlim=xlim,ylim=ylim,truth=truth))
}
