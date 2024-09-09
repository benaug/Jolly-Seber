#This is an SCR version with fixed activity centers, spatial density covariates,
#and an individual survival covariate.

#Recruits are assumed to have same spatial distribution as the alive population
#Survival is assumed to be at random, so spatial distribution of deaths is
#proportional to the spatial distribution of the alive population

#Some notes:
#1) Dimensions probably prevent this from working with only 2 years of data (N.recruit, N.survive are not vectors)
#Need to modify custom updates in this case
#2) Object names that cannot be changed in the nimble model without changes in custom updates:
#N, N.recruit, N.survive, ER, lambda.y1, z.start, z.stop, z.obs, pd
#phi[i,g] (must be of dimension M x n.year),
#Poisson assumptions on N[1] and N.recruit (but can include overdispersion with random effects)
#y can change dimension (e.g., for SCR), but need to account for that in defining "y.nodes"
#below to add custom updates. I think that is all...

library(nimble)
library(coda)
source("sim.JS.SCR.Dcov.R")
source("Nimble Model JS-SCR-Dcov.R")
source("Nimble Functions JS-SCR-Dcov.R") #contains custom distributions and updates
source("sSampler Dcov.R") # activity center sampler that proposes from prior when z.super=0.
source("mask.check.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

n.year <- 4 #number of years
gamma <- rep(0.2,n.year-1) #yearly per-capita recruitment
beta0.phi <- qlogis(0.85) #survival intercept
beta1.phi <- 0.5 #phi response to individual covariate
p0 <- rep(0.1,n.year) #yearly detection probabilities at activity center
sigma <- rep(0.5,n.year) #yearly detection function scale
K <- rep(10,n.year) #yearly sampling occasions

buff <- 3 #state space buffer. Buffers maximal x and y dimensions of X below across years
X <- vector("list",n.year) #one trapping array per year
for(g in 1:n.year){ #using same trapping array every year here
  X[[g]] <- as.matrix(expand.grid(3:11,3:11))
}

### Habitat Covariate stuff###
#buffer maximal trap extent
X.all <- matrix(NA,nrow=0,ncol=2)
for(g in 1:n.year){
  X.all <- rbind(X.all,X[[g]])
}

#get x and y extent by buffering state space
xlim <- range(X.all[,1])+c(-buff,buff)
ylim <- range(X.all[,2])+c(-buff,buff)
#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
x.shift <- xlim[1]
y.shift <- ylim[1]
xlim <- xlim-x.shift
ylim <- ylim-y.shift
for(g in 1:n.year){
  X[[g]][,1] <- X[[g]][,1]-x.shift
  X[[g]][,2] <- X[[g]][,2]-y.shift
}
X.all[,1] <- X.all[,1]-x.shift
X.all[,2] <- X.all[,2]-y.shift

res <- 0.20 #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#create a density covariate
D.cov <- rep(NA,n.cells)
for(c in 1:n.cells){
  D.cov[c] <-  7*dSS[c,1] - 0.5*dSS[c,1]^2 + 7*dSS[c,2] - 0.5*dSS[c,2]^2
}
D.cov <- as.numeric(scale(D.cov))

image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="Covariate Value")
points(X.all,pch=4,cex=0.75,col="lightblue")


#Additionally, maybe we want to exclude "non-habitat"
#just removing the corners here for simplicity
dSS.tmp <- dSS - res/2 #convert back to grid locs
InSS <- rep(1,length(D.cov))
InSS[dSS.tmp[,1]<2&dSS.tmp[,2]<2] <- 0
InSS[dSS.tmp[,1]<2&dSS.tmp[,2]>12] <- 0
InSS[dSS.tmp[,1]>12&dSS.tmp[,2]<2] <- 0
InSS[dSS.tmp[,1]>12&dSS.tmp[,2]>12] <- 0

image(x.vals,y.vals,matrix(InSS,n.cells.x,n.cells.y),main="Habitat")

#Density covariates
D.beta0 <- -2
D.beta1 <- 2
#what is implied expected N in state space?
lambda.cell <- exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

#simulate some data
data <- sim.JS.SCR.Dcov(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InSS=InSS,
            gamma=gamma,n.year=n.year,beta0.phi=beta0.phi,beta1.phi=beta1.phi,
            p0=p0,sigma=sigma,X=X,K=K,xlim=xlim,ylim=ylim,res=res)

#visualize realized activity centers
image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density")
points(X.all,pch=4,cex=0.75)
points(data$truth$s,pch=16)

#function to test for errors in mask set up. 
mask.check(dSS=data$dSS,cells=data$cells,n.cells=data$n.cells,n.cells.x=data$n.cells.x,
           n.cells.y=data$n.cells.y,res=data$res,xlim=data$xlim,ylim=data$ylim,
           x.vals=data$x.vals,y.vals=data$y.vals)


##Initialize##
M <- 300 #data augmentation level. Check N.super posterior to make sure it never hits M
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
s.init <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
idx <- which(rowSums(y.nim)>0) #switch for those actually caught
for(i in idx){
  trps <- matrix(0,nrow=0,ncol=2) #get locations of traps of capture across years for ind i
  for(g in 1:n.year){
    if(sum(y.nim[i,g,])>0){
      trps.g <- matrix(X.nim[g,which(y.nim[i,g,]>0),],ncol=2,byrow=FALSE)
      trps <- rbind(trps,trps.g)
    }
  }
  if(nrow(trps)>1){
    s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
  }else{
    s.init[i,] <- trps
  }
}

#If using a habitat mask, move any s's initialized in non-habitat above to closest habitat
e2dist  <-  function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}
getCell  <-  function(s,res,cells){
  cells[trunc(s[1]/res)+1,trunc(s[2]/res)+1]
}
alldists <- e2dist(s.init,data$dSS)
alldists[,data$InSS==0] <- Inf
for(i in 1:M){
  this.cell <- data$cells[trunc(s.init[i,1]/data$res)+1,trunc(s.init[i,2]/data$res)+1]
  if(data$InSS[this.cell]==0){
    cands <- alldists[i,]
    new.cell <- which(alldists[i,]==min(alldists[i,]))
    s.init[i,] <- data$dSS[new.cell,]
  }
}

#constants for Nimble
#might want to center D.cov here. Simulated D.cov in this testscript is already effectively centered.
constants <- list(n.year=n.year, M=M, J=J, xlim=xlim, ylim=ylim, K1D=K1D,
                  D.cov=data$D.cov,cellArea=data$cellArea,n.cells=data$n.cells,
                  res=data$res)
#inits for Nimble
Niminits <- list(N=N.init,N.survive=N.survive.init,N.recruit=N.recruit.init,
                 ER=N.recruit.init,N.super=N.super.init,z.super=z.super.init,
                 z=z.init,z.start=z.start.init,z.stop=z.stop.init,
                 s=s.init,beta0.phi=0,beta1.phi=0,D0=N.init[1]/(sum(InSS)*res^2),D.beta1=0)

#data for Nimble
dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are
Nimdata <- list(y=y.nim,phi.cov=phi.cov.data,X=X.nim,
                dummy.data=dummy.data,cells=cells,InSS=data$InSS)

# set parameters to monitor
parameters <- c('N','gamma','N.recruit','N.survive','N.super','lambda.y1',
                'beta0.phi','beta1.phi','phi.cov.mu','phi.cov.sd','p0','sigma',
                'D0','D.beta1')
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
config.nodes <- c('beta0.phi','beta1.phi','gamma',paste('phi.cov[',cov.up,']'),
               'phi.cov.mu','phi.cov.sd','p0','sigma','D0','D.beta1')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      nodes=config.nodes,useConjugacy = FALSE)

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

#activity center sampler. This sampler tunes activity centers when z.super[i]=1 and
#draws from the prior otherwise.
# conf$removeSampler(paste0("s[1:",M,", 1:2]")) #dont need to remove if not assigned by nimble above
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSamplerDcov',control=list(i=i,res=data$res,n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
                                                 xlim=data$xlim,ylim=data$ylim),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}


#optional (but recommended!) blocking 
conf$removeSampler(c("beta0.phi"))
conf$removeSampler(c("beta1.phi"))
conf$addSampler(target = c("beta0.phi","beta1.phi"),type = 'RW_block',
                control = list(adaptive=TRUE),silent = TRUE)
conf$addSampler(target = c("D0","D.beta1"),
                type = 'RW_block',control=list(adaptive=TRUE),silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
time1 <- end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2 <- end.time-start.time2 # post-compilation run time

mvSamples <-  as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[-c(1:500),]))

#reminder what the targets are
data$N
data$N.recruit
data$N.survive
data$N[1]+sum(data$N.recruit) #N.super


#Some sanity checks I used during debugging. Just checking that final
#model states match between z and N objects

# #check N count
# N.count <- rep(NA,n.year)
# for(g2 in 1:n.year){
#   N.count[g2] <- sum(Cmodel$z[Cmodel$z.super==1,g2]==1)
# }
# N.count
# Cmodel$N
# 
# #check N.recruit count
# N.count <- rep(NA,n.year-1)
# for(g2 in 2:n.year){
#   N.count[g2-1] <- sum(Cmodel$z[Cmodel$z.super==1,g2-1]==0&Cmodel$z[Cmodel$z.super==1,g2]==1)
# }
# N.count
# Cmodel$N.recruit
# 
# #check N.survive count
# N.count <- rep(NA,n.year-1)
# for(g2 in 2:n.year){
#   N.count[g2-1] <- sum(Cmodel$z[Cmodel$z.super==1,g2-1]==1&Cmodel$z[Cmodel$z.super==1,g2]==1)
# }
# N.count
# Cmodel$N.survive
# 
# #are individual z's consistent with their z.start and z.stop?
# all(apply(Cmodel$z,1,function(x){min(which(x==1))})==Cmodel$z.start)
# all(apply(Cmodel$z,1,function(x){max(which(x==1))})==Cmodel$z.stop)
# 
# #zombie check
# for(i in 1:M){
#   z.tmp <- Cmodel$z[i,]
#   z.on <- which(z.tmp==1)
#   first <- z.on[1]
#   last <- z.on[length(z.on)]
#   if(length(z.on)>1){
#     if(any(z.tmp[first:last]==0)){
#       stop("rawr!")
#     }
#   }
# }