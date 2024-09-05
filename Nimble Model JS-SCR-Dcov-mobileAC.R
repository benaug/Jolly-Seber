NimModel <- nimbleCode({
  #Density covariates
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  # D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  # D.intercept <- exp(D.beta0)*cellArea
  D.intercept <- D0*cellArea
  #First year density model
  for(c in 1:n.cells) {
    lambda.cell[c] <- InSS[c]*exp(D.beta1*D.cov[c]) #separate this component so s's do not depend on D.intercept
    lambda.y1.cell[1,c] <- D.intercept*lambda.cell[c] #expected N in cell c
    pi.cell[c] <- lambda.cell[c] / pi.denom #expected proportion of total N in cell c
  }
  pi.denom <- sum(lambda.cell[1:n.cells])
  
  ##Abundance##
  lambda.y1 <- D.intercept*pi.denom #Expected starting population size
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
  
  #Survival individual covariates
  phi.cov.mu ~ dunif(-10, 10) #phi individual covariate mean prior
  phi.cov.sd ~ T(dt(mu=0, sigma=1, df=7), 0, Inf) #phi individual covariate sd prior
  #RSF coefficient for activity center relocation
  rsf.beta ~ dnorm(0,sd=10)
  sigma.move ~ dunif(0,100) #activity center relocation spatial scale parameter (BVN sd)
  #Resource selection function evaluated across all cells for activity center relocation
  # rsf[1:n.cells] <- exp(rsf.beta*D.cov[1:n.cells]) #if not excluding non-habitat
  rsf[1:n.cells] <- exp(rsf.beta*D.cov[1:n.cells] -5000*OutSS[1:n.cells])
  for(i in 1:M){
    phi.cov[i] ~ dnorm(phi.cov.mu, sd = phi.cov.sd)
    #1st year activity centers
    s[i,1,1] ~ dunif(xlim[1],xlim[2])
    s[i,1,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i,1] <- cells[trunc(s[i,1,1]/res)+1,trunc(s[i,1,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    dummy.data[i] ~ dCell(pi.cell[s.cell[i,1]],InSS[s.cell[i,1]])
    #subsequent year activity center - movement with resource selection
    for(g in 2:n.year){
      #Individual availability distributions conditioned on cells, bivariate Normal centered on activity center
      avail.dist[i,g-1,1:n.cells] <- getAvail(s=s[i,g-1,1:2],sigma=sigma.move,res=res,x.vals=x.vals[1:n.cells.x],
                                          y.vals=y.vals[1:n.cells.y],n.cells.x=n.cells.x,n.cells.y=n.cells.y)
      #Individual use distributions - multiply rsf and available distribution, normalize. not trimmed
      use.dist[i,g-1,1:n.cells] <- getUse(rsf=rsf[1:n.cells],avail=avail.dist[i,g-1,1:n.cells])
      s.cell[i,g] ~ dcat(use.dist[i,g-1,1:n.cells]) #likelihood of this cell being selected
      s.xlim[i,g,1:2] <- dSS[s.cell[i,g],1] + c(-res,res)/2
      s.ylim[i,g,1:2] <- dSS[s.cell[i,g],2] + c(-res,res)/2
      #continuous activity center location likelihood conditioned on the cell
      s[i,g,1] ~ T(dnorm(s[i,g-1,1],sd=sigma.move),s.xlim[i,g,1],s.xlim[i,g,2])
      s[i,g,2] ~ T(dnorm(s[i,g-1,2],sd=sigma.move),s.ylim[i,g,1],s.ylim[i,g,2])
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
#4) activity centers. Propose entire individual movement trajectory from prior when z.super=0