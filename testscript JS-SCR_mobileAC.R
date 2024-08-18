#This is an SCR version with mobile activity centers, and an individual survival covariate.

library(nimble)
library(coda)
library(truncnorm) #for data simulator
source("sim.JS.SCR.R")
source("Nimble Model JS-SCR_mobileAC.R")
source("Nimble Functions JS-SCR_mobileAC.R") #contains custom distributions and updates
source("sSampler Multi.R") # activity center sampler that proposes from prior when z.super=0.
#this one works for fixed activity centers over years only

n.year <- 4 #number of years
lambda.y1 <- 100 #expected N in year 1
gamma <- rep(0.2,n.year-1) #yearly per-capita recruitment
beta0.phi <- qlogis(0.85) #survival intercept
beta1.phi <- 0.5 #phi response to individual covariate
p0 <- rep(0.1,n.year) #yearly detection probabilities at activity center
sigma <- rep(0.5,n.year) #yearly detection function scale
sigma.move <- 0.75 #movement sigma, fixed over primary periods
K <- rep(10,n.year) #yearly sampling occasions

buff <- 2 #state space buffer. Buffers maximal x and y dimensions of X below across years
X <- vector("list",n.year) #one trapping array per year
for(g in 1:n.year){ #using same trapping array every year here
  X[[g]] <- as.matrix(expand.grid(3:11,3:11))
}

data <- sim.JS.SCR(lambda.y1=lambda.y1,gamma=gamma,n.year=n.year,
            beta0.phi=beta0.phi,beta1.phi=beta1.phi,
            p0=p0,sigma=sigma,X=X,buff=buff,K=K,sigma.move=sigma.move)


##Initialize##
M <- 200 #data augmentation level. Check N.super posterior to make sure it never hits M
N.super.init <- nrow(data$y)
X <- data$X #pull X from data (won't be in environment if not simulated directly above)
if(N.super.init > M) stop("Must augment more than number of individuals captured")
J <- unlist(lapply(X,nrow)) #traps per year
J.max <- max(J)

y.nim <- array(0,dim=c(M,n.year,J.max))
y.nim[1:N.super.init,1:n.year,1:J.max] <- data$y #all these guys must be observed
#initialize z, start with observed guys
z.init <- 1*(y.nim>0)
z.init <- matrix(0,M,n.year)
z.start.init <- z.stop.init <- rep(NA,M)
y.nim2D <- apply(y.nim,c(1,2),sum)
for(i in 1:N.super.init){
  det.idx <- which(y.nim2D[i,]>0)
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
z.obs <- 1*(rowSums(y.nim)>0) #indicator for "ever observed"

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

#remaining SCR stuff to initialize
#put X in ragged array
#also make K1D, year by trap operation history, as ragged array.
X.nim <- array(0,dim=c(n.year,J.max,2))
K1D <- matrix(0,n.year,J.max)
for(g in 1:n.year){
  X.nim[g,1:J[g],1:2] <- X[[g]]
  K1D[g,1:J[g]] <- rep(K[g],J[g])
}

#pull out state space with buffer around maximal trap dimensions
xlim <- data$xlim
ylim <- data$ylim
s.init <- array(NA,dim=c(M,n.year,2))
for(g in 1:n.year){
  s.init[,g,]=cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
  idx=which(rowSums(y.nim[,g,])>0) #switch for those actually caught
  for(i in idx){
    trps <- matrix(X.nim[g,which(y.nim[i,g,]>0),],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s.init[i,g,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s.init[i,g,] <- trps
    }
  }
}

#constants for Nimble
constants <- list(n.year=n.year, M=M, J=J, xlim=xlim, ylim=ylim, K1D=K1D)
#inits for Nimble
Niminits <- list(N=N.init,N.survive=N.survive.init,N.recruit=N.recruit.init,
                 z=z.init,z.start=z.start.init,z.stop=z.stop.init,
                 ER=N.recruit.init,N.super=N.super.init,z.super=z.super.init,
                 s=s.init,beta0.phi=0,beta1.phi=0)

#data for Nimble
Nimdata <- list(y=y.nim,phi.cov=phi.cov.data,X=X.nim)

# set parameters to monitor
parameters <- c('N','gamma','N.recruit','N.survive','N.super',
                'lambda.y1','beta0.phi','beta1.phi',
                'phi.cov.mu','phi.cov.sd','p0','sigma','sigma.move')
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
config.nodes <- c('beta0.phi','beta1.phi','gamma','lambda.y1',paste('phi.cov[',cov.up,']'),
               'phi.cov.mu','phi.cov.sd','p0','sigma','sigma.move')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      nodes=config.nodes,useConjugacy = TRUE)

###*required* sampler replacements###
z.super.ups <- round(M*0.2) #how many z.super update proposals per iteration? 
#20% of M seems reasonable, but optimal will depend on data set
#loop here bc potentially different numbers of traps to vectorize in each year
y.nodes <- pd.nodes <- c()
for(g in 1:n.year){
  y.nodes <- c(y.nodes,Rmodel$expandNodeNames(paste0("y[1:",M,",",g,",1:",J[g],"]"))) #if you change y structure, change here
  pd.nodes <- c(pd.nodes,Rmodel$expandNodeNames(paste0("pd[1:",M,",",g,",1:",J[g],"]"))) #if you change y structure, change here
}
N.nodes <- Rmodel$expandNodeNames(paste0("N"))
N.survive.nodes <- Rmodel$expandNodeNames(paste0("N.survive[1:",n.year-1,"]"))
N.recruit.nodes <- Rmodel$expandNodeNames(paste0("N.recruit[1:",n.year-1,"]"))
ER.nodes <- Rmodel$expandNodeNames(paste0("ER[1:",n.year-1,"]"))
z.nodes <- Rmodel$expandNodeNames(paste0("z[1:",M,",1]"))
calcNodes <- c(N.nodes,N.recruit.nodes,y.nodes,z.nodes) #the ones that need likelihoods updated in mvSaved
conf$addSampler(target = c("z"),
                type = 'zSampler',control = list(M=M,n.year=n.year,J=J,
                                                 z.obs=z.obs,z.super.ups=z.super.ups,
                                                 y.nodes=y.nodes,pd.nodes=pd.nodes,N.nodes=N.nodes,
                                                 z.nodes=z.nodes,ER.nodes=ER.nodes,
                                                 N.survive.nodes=N.survive.nodes,
                                                 N.recruit.nodes=N.recruit.nodes,
                                                 y2D=y.nim2D,calcNodes=calcNodes), silent = TRUE)

#activity center sampler. There are 3 samplers here for these cases
#1) z.super=1 and z=1, sSampler1 uses Metropolis-Hastings
#2) z.super=1 and z=0, sSampler2 uses Metropolis-Hastings with proposal sd tuned separately from above
#3) z.super=0, sSampler3 simulates entire new trajectory from priors.
for(i in 1:M){
  for(g in 1:n.year){
    conf$addSampler(target = paste0("s[",i,",",g,",1:2]"),
                    type = 'sSampler1',control=list(i=i,g=g,xlim=xlim,ylim=ylim,scale=1),silent = TRUE)
    conf$addSampler(target = paste0("s[",i,",",g,",1:2]"),
                    type = 'sSampler2',control=list(i=i,g=g,xlim=xlim,ylim=ylim,scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
  conf$addSampler(target = paste0("s[",i,",1:",n.year,",1:2]"),
                  type = 'sSampler3',control=list(i=i,xlim=xlim,ylim=ylim,n.year=n.year,scale=1),silent = TRUE)
}


#optional (but recommended!) blocking 
#do i need this in this model? probably
conf$removeSampler(c("beta0.phi"))
conf$removeSampler(c("beta1.phi"))
conf$addSampler(target = c("beta0.phi","beta1.phi"),type = 'RW_block',
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
data$N[1]+sum(data$N.recruit) #N.super


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
N.count <- rep(NA,n.year-1)
for(g2 in 2:n.year){
  N.count[g2-1] <- sum(Cmodel$z[Cmodel$z.super==1,g2-1]==0&Cmodel$z[Cmodel$z.super==1,g2]==1) 
}
N.count
Cmodel$N.recruit

#check N.survive count
N.count <- rep(NA,n.year-1)
for(g2 in 2:n.year){
  N.count[g2-1] <- sum(Cmodel$z[Cmodel$z.super==1,g2-1]==1&Cmodel$z[Cmodel$z.super==1,g2]==1) 
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