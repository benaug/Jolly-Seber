NimModel <- nimbleCode({
  ##Abundance##
  lambda.y1 ~ dunif(0,1000) #Expected starting population size
  N[1] ~ dpois(lambda.y1) #Realized starting population size
  for(g in 2:n.year){
    N[g] <- N.survive[g-1] + N.recruit[g-1] #yearly abundance
    #N.recruit and N.survive information also contained in z/z.start + z.stop
    #N.recruit has distributions assigned below, but survival distributions defined on z
  }
  N.super <- N[1] + sum(N.recruit[1:(n.year-1)]) #size of superpopulation
  
  #Recruitment
  for(g in 1:(n.year-1)){
    gamma[g] ~ dunif(0,2) # yearly recruitment priors
    ER[g] <- N[g]*gamma[g] #yearly expected recruits
    N.recruit[g] ~ dpois(ER[g]) #yearly realized recruits
  }
  
  #Individual covariates
  phi.cov.mu ~ dunif(-10, 10) #phi individual covariate mean prior
  phi.cov.sd ~ T(dt(mu=0, sigma=1, df=7), 0, Inf) #phi individual covariate sd prior
  sigma.move ~ dunif(0,2)
  for(i in 1:M){
    phi.cov[i] ~ dnorm(phi.cov.mu, sd = phi.cov.sd)
    #1st year ACs
    s[i,1,1] ~ dunif(xlim[1],xlim[2])
    s[i,1,2] ~ dunif(ylim[1],ylim[2])
    for(g in 2:n.year){
      #If you change movement distribution, need to modify custom s updates (3rd one for z.super=0)
      s[i,g,1] ~ T(dnorm(s[i,g-1,1],sd=sigma.move),xlim[1],xlim[2])
      s[i,g,2] ~ T(dnorm(s[i,g-1,2],sd=sigma.move),ylim[1],ylim[2])
    }
  }
  
  #Survival (phi must have M x n.year - 1 dimension for custom updates to work)
  #without individual or year effects, use for loop to plug into phi[i,g]
  beta0.phi ~ dlogis(0,1)
  beta1.phi ~ dnorm(0, sd=10) #individual covariate effect on survival
  for(i in 1:M){
    for(g in 1:(n.year-1)){#plugging same individual phi's into each year for custom update
      logit(phi[i,g]) <- beta0.phi + beta1.phi*phi.cov[i] #individual by year survival
    }
    #survival likelihood (bernoulli) that only sums from z.start to z.stop
    z[i,1:n.year] ~ dSurvival(phi=phi[i,1:(n.year-1)],z.start=z.start[i],z.stop=z.stop[i])
  }
  
  ##Detection##
  sigma ~ dunif(0,10) #fixing sigma across years
  for(g in 1:n.year){
    p0[g] ~ dunif(0,1) #p0 varies by year
    for(i in 1:M){ #only compute pd and y when z.super[i]=1&z[i,g]=1
      pd[i,g,1:J[g]] <- GetDetectionProb(s = s[i,g,1:2], X = X[g,1:J[g],1:2],
                                         J=J[g], sigma=sigma, p0=p0[g], z=z[i,g], z.super=z.super[i])
      y[i,g,1:J[g]] ~ dBinomialVector(pd = pd[i,g,1:J[g]], K = K1D[g,1:J[g]],
                                      z = z[i,g], z.super=z.super[i]) #vectorized obs mod
    }
  }
})

#custom updates:
#1) for detected individuals: update z.start, then update z.stop
#2) for undetected individuals: update entire z vectors
#3) N.super/z.super update