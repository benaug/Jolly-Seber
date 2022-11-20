sim.JS <- function(lambda.y1=NA,gamma=NA,beta0.phi=NA,beta1.phi=NA,beta2.phi=NA,
                   p=NA,n.year=NA,K=NA){
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
    phi[,g-1]=plogis(beta0.phi+beta1.phi[g-1]+cov*beta2.phi)
    
    idx=which(z[,g-1]==1)
    z[idx,g]=rbinom(length(idx),1,phi[idx,g-1])
    N.survive[g]=sum(z[,g-1]==1&z[,g]==1)
    N[g]=N.recruit[g]+N.survive[g]
  }
  
  if(any(N.recruit[2:n.year]+N.survive[2:n.year]!=N[2:n.year]))stop("Simulation bug")
  if(any(colSums(z)!=N))stop("Simulation bug")
  
  #plot to see if sim values realistic
  hist(phi,main="dist of phi across ind and time")
  
  birth.year=apply(z,1,function(x){match(1,x)})
  death.year=10-apply(z,1,function(x){match(1,rev(x))})+1
  lifetimes=death.year-birth.year
  
  #detection
  y=z*0
  for(g in 1:n.year){
      y[,g]=rbinom(nrow(z),K[g],p[g]*z[,g])
  }
  
  #store true data for model buildling/debugging
  truth=list(y=y,cov=cov,N=N,N.recruit=N.recruit,N.survive=N.survive,z=z)
  
  #discard undetected individuals
  keep.idx=which(rowSums(y)>0)
  y=y[keep.idx,]
  cov=cov[keep.idx]
  return(list(y=y,cov=cov,N=N,N.recruit=N.recruit,N.survive=N.survive,truth=truth))
}