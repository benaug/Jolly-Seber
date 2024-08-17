NimModel <- nimbleCode({
  ##Abundance##
  lambda.y1.M ~ dunif(0,1000) #Expected starting population size, males
  lambda.y1.F ~ dunif(0,1000) #Expected starting population size, females
  N.M[1] ~ dpois(lambda.y1.M) #Realized starting population size
  N.F[1] ~ dpois(lambda.y1.F) #Realized starting population size
  N[1] <- N.M[1] + N.F[1]
  for(g in 2:n.year){
    N[g] <- N.survive[g-1] + N.recruit[g-1] #yearly total abundance
    N.M[g] <- N.survive.M[g-1] + N.recruit.M[g-1] #yearly male abundance
    N.F[g] <- N.survive.F[g-1] + N.recruit.F[g-1] #yearly female abundance
    #sex-specific N.recruit and N.survive information also contained in z, z.start, z.stop, and sex
    #sex-specific N.recruit has distributions assigned below, but survival distributions defined on z
  }
  N.super <- N[1] + sum(N.recruit[1:(n.year-1)]) #size of superpopulation
  
  #Sex-specific Recruitment
  gamma.male ~ dunif(0,2)
  gamma.female ~ dunif(0,2)
  for(g in 1:(n.year-1)){
    # gamma.sex[1,g] ~ dunif(0,2) # yearly male recruitment priors
    # gamma.sex[2,g] ~ dunif(0,2) # yearly female recruitment priors
    gamma.sex[1,g] <- gamma.male
    gamma.sex[2,g] <- gamma.female
    ER.M[g] <- N[g]*gamma.sex[1,g] #yearly male expected recruits per total N
    ER.F[g] <- N[g]*gamma.sex[2,g] #yearly female expected recruits per total N
    N.recruit.M[g] ~ dpois(ER.M[g]) #yearly male realized recruits
    N.recruit.F[g] ~ dpois(ER.F[g]) #yearly female realized recruits
  }
  
  #Individual covariates. individual sex info derived from recruitment model
  for(i in 1:M){
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
  }

  #Survival (phi must have M x n.year - 1 dimension for custom updates to work)
  #without individual or year effects, use for loop to plug into phi[i,g]
  phi.sex[1] ~ dunif(0,1)
  phi.sex[2] ~ dunif(0,1)
  
  for(i in 1:M){
    for(g in 1:(n.year-1)){#plugging same individual phi's into each year (phi: M x n.year-1 expected by custom update)
      phi[i,g] <- phi.sex[sex[i]+1]
    }
    #survival likelihood (bernoulli) that only sums from z.start to z.stop
    z[i,1:n.year] ~ dSurvival(phi=phi[i,1:(n.year-1)],z.start=z.start[i],z.stop=z.stop[i])
  }

  ##Detection##
  #fixing sex-specific sigma across years
  sigma.sex[1] ~ dunif(0,10) #male sigma
  sigma.sex[2] ~ dunif(0,10) #female sigma
  for(g in 1:n.year){
    #sex-specific p0 varies by year
    p0.sex[1,g] ~ dunif(0,1) #male p0
    p0.sex[2,g] ~ dunif(0,1) #female p0
    for(i in 1:M){ #only compute pd and y likelihood when z.super[i]=1&z[i,g]=1
      pd[i,g,1:J[g]] <- GetDetectionProb(s = s[i,1:2], X = X[g,1:J[g],1:2],
                                         J=J[g], sigma=sigma.sex[sex[i]+1], p0=p0.sex[sex[i]+1,g], z=z[i,g], z.super=z.super[i])
      y[i,g,1:J[g]] ~ dBinomialVector(pd = pd[i,g,1:J[g]], K = K1D[g,1:J[g]],
                                      z = z[i,g], z.super=z.super[i]) #vectorized obs mod
    }
  }
})

#custom updates:
#1) for detected individuals: update z.start, then update z.stop
#2) for undetected individuals: update entire z vectors
#3) N.super/z.super update
#4) sex update: update sex-specific N structures with sex[i] update