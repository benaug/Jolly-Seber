sim.JS <- function(lambda.y1=NA,gamma=NA,beta0.phi=NA,beta1.phi=NA,beta2.phi=NA,
                   p=NA,n.year=NA,K=NA){
  #Population dynamics
  N=rep(NA,n.year)
  N.recruit=N.survive=ER=rep(NA,n.year-1)
  N[1]=rpois(1,lambda.y1)
  
  #Easiest to increase dimension of z as we simulate bc size not known in advance.
  z=matrix(0,N[1],n.year)
  z[1:N[1],1]=1
  cov=rnorm(N[1],0,1) #simulate ind survival covariate for 1st year guys
  phi=matrix(NA,N[1],n.year-1)
  for(g in 2:n.year){
    #Simulate recruits
    ER[g-1]=N[g-1]*gamma[g-1]
    N.recruit[g-1]=rpois(1,ER[g-1])
    #add recruits to z
    z.dim.old=length(cov)
    z=rbind(z,matrix(0,nrow=N.recruit[g-1],ncol=n.year))
    z[(z.dim.old+1):(z.dim.old+N.recruit[g-1]),g]=1
    cov=c(cov,rep(NA,N.recruit[g-1]))
    cov[(z.dim.old+1):(z.dim.old+N.recruit[g-1])]=rnorm(N.recruit[g-1],0,1) #simulate survival cov values for new recruits
    
    #Simulate survival
    phi=rbind(phi,matrix(NA,nrow=N.recruit[g-1],ncol=n.year-1))
    phi[,g-1]=plogis(beta0.phi+cov*beta1.phi)
    
    idx=which(z[,g-1]==1)
    z[idx,g]=rbinom(length(idx),1,phi[idx,g-1])
    N.survive[g-1]=sum(z[,g-1]==1&z[,g]==1)
    N[g]=N.recruit[g-1]+N.survive[g-1]
  }
  
  if(any(N.recruit+N.survive!=N[2:n.year]))stop("Simulation bug")
  if(any(colSums(z)!=N))stop("Simulation bug")
  
  #plot to see if sim values realistic
  hist(phi,main="Distribution of Individual Phi")
  
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