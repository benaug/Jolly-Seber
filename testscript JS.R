#Some notes:
#1) Dimensions probably prevent this from working with only 2 years of data (N.recruit, N.survive are not vectors)
#Need to modify custom updates in this case
#2) Object names that cannot be changed in the nimble model without changes in custom updates:
#N, N.recruit, N.survive, ER, lambda.y1, z.start, z.stop, z.obs, 
#phi[i,g] (must be of dimension M x n.year),
#Poisson assumptions on N[1] and N.recruit (but can include overdispersion with random effects)
#y can change dimension (e.g., for SCR), but need to account for that in defining "y.nodes"
#below to add custom updates. I think that is all...

library(nimble)
library(coda)
source("sim.JS.R")
source("Nimble Model JS.R")
source("Nimble Functions JS.R") #contains custom distributions and updates

#make sure to run this line or the MCMC sampler will not work!
nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)
nimbleOptions('MCMCjointlySamplePredictiveBranches') 

n.year <- 5 #number of years
lambda.y1 <- 200 #expected N in year 1
gamma<- rep(0.2,n.year) #yearly per-capita recruitment (there are only n.year-1 recruitment parameters, data simulator ignores 1st entry, need to fix this.)
beta0.phi <- qlogis(0.85) #survival intercept
#simulating yearly survival offsets here, but really need large populations or simpler model to estimate with much
#precision, at least with effect sizes simulated here (can increase the simulation variance on next line)
beta1.phi <- c(0,rnorm(n.year-2,0,0.2)) #survival year offsets, year 1 is intercept. must be of length n.year-1
plogis(beta0.phi+beta1.phi) #implied yearly survival at cov value of 0
beta2.phi <- 0.5 #phi response to individual covariate
p <- rep(0.1,n.year) #yearly detection probabilities
K <- rep(10,n.year) #yearly sampling occasions

data <- sim.JS(lambda.y1=lambda.y1,gamma=gamma,
            beta0.phi=beta0.phi,beta1.phi=beta1.phi,
            beta2.phi=beta2.phi,p=p,n.year=n.year,K=K)


##Initialize##
M <- 500 #data augmentation level. Check N.super posterior to make sure it never hits M
N.super.init <- nrow(data$y)
if(N.super.init > M) stop("Must augment more than number of individuals captured")
y.init <- matrix(0,M,n.year)
y.init[1:N.super.init,] <- data$y #all these guys must be observed
#initialize z, start with observed guys
z.init <- 1*(y.init>0)
z.init <- matrix(0,M,n.year)
z.start.init <- z.stop.init <- rep(NA,M)
for(i in 1:N.super.init){
  det.idx <- which(y.init[i,]>0)
  if(length(det.idx)>0){ #some all 0 histories if init from truth
    z.start.init[i] <- min(det.idx)
    z.stop.init[i] <- max(det.idx)
    z.init[i,z.start.init[i]:z.stop.init[i]] <- 1
  }
}
#now initialize augmented z's
for(i in (N.super.init+1):M){
  start <- sample(1:n.year,1) #random recruit year
  if(start<(n.year-1)){ #random death year
    stop <- sample(start:n.year,1)
  }else{ #unless recruited one year before end
    stop <- n.year
  }
  z.init[i,start:stop] <- 1
  z.start.init[i] <- start
  z.stop.init[i] <- stop
}
z.super.init <- c(rep(1,N.super.init),rep(0,M-N.super.init))
z.obs.init <- 1*(rowSums(y.init)>0) #indicator for "ever observed"

#initialize N structures from z.init
N.init <- colSums(z.init[z.super.init==1,])
N.survive.init <- N.recruit.init <- rep(NA,n.year-1)
for(g in 2:n.year){
  N.survive.init[g-1] <- sum(z.init[,g-1]==1&z.init[,g]==1&z.super.init==1)
  N.recruit.init[g-1] <- N.init[g]-N.survive.init[g-1]
}

#now individual covariate
phi.cov.data <- c(data$cov,rep(NA,M-length(data$cov)))
cov.up <- which(is.na(phi.cov.data)) #which individuals have missing cov values, used below to help nimble assign samplers
#constants for Nimble
constants <- list(n.year=n.year, K=K, M=M)
#inits for Nimble
Niminits <- list(N=N.init,N.survive=N.survive.init,N.recruit=N.recruit.init,
                 ER=N.recruit.init,N.super=N.super.init,z.super=z.super.init,
                 beta0.phi=beta0.phi,beta1.phi=beta1.phi[-1],lambda.y1=lambda.y1) #for demonstration. Don't start at truth in practice.

#data for Nimble
Nimdata <- list(y=y.init,z=z.init,z.start=z.start.init,z.stop=z.stop.init,phi.cov=phi.cov.data)

# set parameters to monitor
parameters <- c('N','gamma','N.recruit','N.survive','p','beta0.phi',
              'beta1.phi','beta2.phi','lambda.y1','N.super','phi.cov.mu','phi.cov.sd')
nt <- 1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)

#OK! what are we doing here? If you just let nimble configure as normal, it will assign incorrect samplers
#to z and N objects. We could then remove them and replace them, but it is much faster to not let nimble
#make the assignments in the first place. So! put all terms with priors in config.nodes here except for
#N.recruit. If you change the model parameters, you will need to make the same changes here. Finally, 
#we have to tell nimble which nodes to assign samplers for for the individual covariate when manually
#instructing nimble which samplers to assign.
config.nodes <- c('beta0.phi','beta1.phi','beta2.phi','gamma','p','lambda.y1',paste('phi.cov[',cov.up,']'),
               'phi.cov.mu','phi.cov.sd')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      nodes=config.nodes,useConjugacy = TRUE)

###*required* sampler replacement###
z.super.ups <- round(M*0.2) #how many z.super update proposals per iteration? 
#20% of M seems reasonable, but optimal will depend on data set
y.nodes <- Rmodel$expandNodeNames(paste0("y[1:",M,",1]")) #if you change y structure, change here
N.nodes <- Rmodel$expandNodeNames(paste0("N"))
N.survive.nodes <- Rmodel$expandNodeNames(paste0("N.survive[1:",n.year-1,"]"))
N.recruit.nodes <- Rmodel$expandNodeNames(paste0("N.recruit[1:",n.year-1,"]"))
ER.nodes <- Rmodel$expandNodeNames(paste0("ER[1:",n.year-1,"]"))
z.nodes <- Rmodel$expandNodeNames(paste0("z[1:",M,",1]"))
calcNodes <- c(N.nodes,N.recruit.nodes,y.nodes,z.nodes) #the ones that need likelihoods updated in mvSaved
conf$addSampler(target = c("z"),
                type = 'zSampler',control = list(M=M,n.year=n.year,z.obs=z.obs.init,z.super.ups=z.super.ups,
                                                 y.nodes=y.nodes,N.nodes=N.nodes,z.nodes=z.nodes,
                                                 ER.nodes=ER.nodes,N.survive.nodes=N.survive.nodes,
                                                 N.recruit.nodes=N.recruit.nodes,
                                                 calcNodes=calcNodes), silent = TRUE)

#optional (but recommended!) blocking 
#strong posterior correlation between intercept and yearly offsets
conf$removeSampler(c("beta0.phi"))
conf$removeSampler(c("beta1.phi"))
conf$addSampler(target = c("beta0.phi","beta1.phi"),type = 'RW_block',
                control = list(adaptive=TRUE),silent = TRUE)
#blocking p should increase efficiency because with the y node vectorized over years, you must evaluate
#all year's detections to update 1 year's p. So blocking replaces n.year likelihood evaluations with 1.
conf$removeSampler(c("p"))
conf$addSampler(target = c("p"),type = 'RW_block',
                control = list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(5000,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
time1 <- end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2 <- end.time-start.time2 # post-compilation run time

mvSamples <-  as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[-c(1:1000),]))

#reminder what the targets are
data$N
data$N.recruit
data$N.survive
data$N[1]+sum(data$N.recruit[-1]) #N.super




#Some sanity checks I used during debugging. Just checking that final
#model states match between z and N objects

#check N count
N.count <- rep(NA,n.year)
for(g2 in 1:n.year){
  N.count[g2] <- sum(Cmodel$z[Cmodel$z.super==1,g2]==1) 
}
N.count
Cmodel$N

#check N.recruit count
N.count <- rep(NA,n.year)
for(g2 in 1:n.year){
  N.count[g2] <- sum(Cmodel$z[Cmodel$z.super==1,g2-1]==0&Cmodel$z[Cmodel$z.super==1,g2]==1) 
}
N.count
Cmodel$N.recruit

#check N.survive count
N.count <- rep(NA,n.year)
for(g2 in 1:n.year){
  N.count[g2] <- sum(Cmodel$z[Cmodel$z.super==1,g2-1]==1&Cmodel$z[Cmodel$z.super==1,g2]==1) 
}
N.count
Cmodel$N.survive

#are individual z's consistent with their z.start and z.stop?
all(apply(Cmodel$z,1,function(x){min(which(x==1))})==Cmodel$z.start)
all(apply(Cmodel$z,1,function(x){max(which(x==1))})==Cmodel$z.stop)


#zombie check
for(i in 1:M){
  z.tmp <- Cmodel$z[i,]
  z.on <- which(z.tmp==1)
  first <- z.on[1]
  last <- z.on[length(z.on)]
  if(length(z.on)>1){
    if(any(z.tmp[first:last]==0)){
      stop("rawr!")
    }
  }
}