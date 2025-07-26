#This is an SCR version with spatial density covariates, RSF activity center movement,
#and an individual survival covariate. Activity centers are in continuous space.
#Activity center relocation uses a proper resource selection model with an availability 
#distribution that is bivariate normal. Could use other availability distributions, bivariate t, etc.

#This version has not been tested. Runs pretty slowly. May be able to improve efficiency more, will look into that.

library(nimble)
library(coda)
source("sim.JS.SCR.Dcov.mobileAC.R")
source("Nimble Model JS-SCR-Dcov-mobileAC.R")
source("Nimble Functions JS-SCR-Dcov-mobileAC.R") #contains custom distributions and updates
source("sSampler Dcov MobileAC.R") # required activity center samplers
source("mask.check.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

n.year <- 5 #number of years
gamma <- rep(0.2,n.year-1) #yearly per-capita recruitment
beta0.phi <- qlogis(0.85) #survival intercept
beta1.phi <- 0.5 #phi response to individual covariate
p0 <- rep(0.1,n.year) #yearly detection probabilities at activity center
sigma <- rep(0.75,n.year) #yearly detection function scale
sigma.move <- 2 #movement sigma, fixed over primary periods
rsf.beta <- 1 #selection coefficient for activity center relocation btwn primary periods
K <- rep(10,n.year) #yearly sampling occasions

buff <- 6 #state space buffer. Buffers maximal x and y dimensions of X below across years
X <- vector("list",n.year) #one trapping array per year
for(g in 1:n.year){ #using same trapping array every year here
  # X[[g]] <- as.matrix(expand.grid(3:11,3:11))
  X[[g]] <- as.matrix(expand.grid(seq(from=3*sigma[1],by=2*sigma[1],length.out=9),
                                  seq(from=3*sigma[1],by=2*sigma[1],length.out=9)))
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

res <- 0.50 #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")

#create a density covariate
#simulate a D.cov, higher cov.pars for large scale cov
set.seed(1320570)
library(geoR)
D.cov <- grf(n.cells,grid=dSS,cov.pars=c(25,25),messages=FALSE)[[2]]
D.cov <- as.numeric(scale(D.cov)) #scale
par(mfrow=c(1,1),ask=FALSE)
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)

image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="Covariate Value")
points(X.all,pch=4,cex=0.75,col="lightblue",lwd=2)


#Additionally, maybe we want to exclude "non-habitat"
#in first year D model
#just removing the corners here for simplicity
dSS.tmp <- dSS - res/2 #convert back to grid locs
InSS <- rep(1,length(D.cov))
InSS[dSS.tmp[,1]<(xlim[1]+1)&dSS.tmp[,2]<(ylim[1]+1)] <- 0
InSS[dSS.tmp[,1]<(xlim[1]+1)&dSS.tmp[,2]>=(ylim[2]-1)] <- 0
InSS[dSS.tmp[,1]>=(xlim[2]-1)&dSS.tmp[,2]<(ylim[1]+1)] <- 0
InSS[dSS.tmp[,1]>=(xlim[2]-1)&dSS.tmp[,2]>=(ylim[2]-1)] <- 0

image(x.vals,y.vals,matrix(InSS,n.cells.x,n.cells.y),main="Habitat")

#Density covariates
D.beta0 <- -2.5
D.beta1 <- 1
#what is implied expected N in state space?
lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

#simulate some data
data <- sim.JS.SCR.Dcov.mobileAC(D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InSS=InSS,
            gamma=gamma,n.year=n.year,beta0.phi=beta0.phi,beta1.phi=beta1.phi,
            p0=p0,sigma=sigma,sigma.move=sigma.move,rsf.beta=rsf.beta,
            X=X,K=K,xlim=xlim,ylim=ylim,res=res)

#visualize realized activity centers in a given year
#compare first and last year to see if/how spatial distribution of activity centers changed over time.
par(mfrow=c(1,1),ask=FALSE)
plot.year <- 1
image(x.vals,y.vals,matrix(data$truth$pi.cell[plot.year,],n.cells.x,n.cells.y),
      main=paste("Expected proportion of realized N in each cell, year", plot.year))
points(X.all,pch=4,cex=0.75)
points(data$truth$s[data$truth$z[,plot.year]==1,plot.year,],pch=16)

# can look at individual by year availability and use distributions
# i <- 1
# g <- 1
# par(mfrow=c(3,1))
# image(x.vals,y.vals,matrix(data$truth$avail.dist[i,g,],n.cells.x,n.cells.y),main="Availability Distribution")
# points(data$truth$s[i,g,1],data$truth$s[i,g,2],pch=16,col="lightblue")
# image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="RSF Cov")
# points(data$truth$s[i,g,1],data$truth$s[i,g,2],pch=16,col="lightblue")
# image(x.vals,y.vals,matrix(data$truth$use.dist[i,g,],n.cells.x,n.cells.y),main="Use Distribution")
# points(data$truth$s[i,g,1],data$truth$s[i,g,2],pch=16,col="lightblue")
# par(mfrow=c(1,1))

#function to test for errors in mask set up. 
mask.check(dSS=data$dSS,cells=data$cells,n.cells=data$n.cells,n.cells.x=data$n.cells.x,
           n.cells.y=data$n.cells.y,res=data$res,xlim=data$xlim,ylim=data$ylim,
           x.vals=data$x.vals,y.vals=data$y.vals)


##Initialize##
M <- 225 #data augmentation level. Check N.super posterior to make sure it never hits M
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

#initialize activity centers. Could be improved, but this should be decent
#1) initialize captured guy ACs in years captured, set to mean of captured years in uncaptured years
s.init <- array(NA,dim=c(M,n.year,2))
for(i in 1:N.super.init){
  cap.years <- which(rowSums(y.nim[i,,])>0) 
  nocap.years <- setdiff(1:n.year,cap.years)
  for(g in cap.years){
    #add small increment so we don't initialize exactly on a cell boundary, can lead to -Inf starting logprob due to BNV inside cell
    trps <- matrix(X.nim[g,which(y.nim[i,g,]>0),],ncol=2,byrow=FALSE) + 0.000001
    if(nrow(trps)>1){
      s.init[i,g,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s.init[i,g,] <- trps
    }
  }
  if(length(nocap.years)>0){
    if(length(cap.years)>1){
      s.tmp <- colMeans(s.init[i,cap.years,])
    }else{
      s.tmp <- s.init[i,cap.years,]
    }
    s.init[i,nocap.years,1] <- s.tmp[1]
    s.init[i,nocap.years,2] <- s.tmp[2]
  }
}

#2) get a starting sigma.move that will have a finite logprob
#use to simulating s.inits for augmented individuals
sigma.move.init <- 0
for(g in 2:n.year){
  tmp.dists <- sqrt((s.init[1:N.super.init,g,1]-s.init[1:N.super.init,g-1,1])^2+(s.init[1:N.super.init,g,2]-s.init[1:N.super.init,g-1,2])^2)
  d.max <- max(c((s.init[1:N.super.init,g,1]-s.init[1:N.super.init,g-1,1])^2,(s.init[1:N.super.init,g,2]-s.init[1:N.super.init,g-1,2])^2))
  sigma.move.tmp <- sqrt(d.max)
  if(sigma.move.tmp>sigma.move.init){
    sigma.move.init <- sigma.move.tmp
  }
}

#3) simulate trajectories for the augmented individuals using sigma.move.init
library(truncnorm)
for(i in (N.super.init+1):M){
  s.init[i,1,] <- cbind(runif(1,xlim[1],xlim[2]), runif(1,ylim[1],ylim[2])) #assign random locations
  for(g in 2:n.year){
    s.init[i,g,1] <- rtruncnorm(1,xlim[1],xlim[2],mean=s.init[i,g-1,1],sd=sigma.move.init)
    s.init[i,g,2] <- rtruncnorm(1,ylim[1],ylim[2],mean=s.init[i,g-1,2],sd=sigma.move.init)
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
for(g in 1:n.year){
  alldists <- e2dist(s.init[,g,],data$dSS)
  alldists[,data$InSS==0] <- Inf
  for(i in 1:M){
    this.cell <- data$cells[trunc(s.init[i,g,1]/data$res)+1,trunc(s.init[i,g,2]/data$res)+1]
    if(data$InSS[this.cell]==0){
      cands <- alldists[i,]
      new.cell <- which(alldists[i,]==min(alldists[i,]))
      s.init[i,g,] <- data$dSS[new.cell,]
    }
  }
}

#get s.cell.inits
s.cell.init <- matrix(NA,M,n.year)
for(i in 1:M){
  for(g in 1:n.year){
    s.cell.init[i,g] <- data$cells[trunc(s.init[i,g,1]/data$res)+1,trunc(s.init[i,g,2]/data$res)+1]
  }
}

#constants for Nimble
#might want to center D.cov here. Simulated D.cov in this testscript is already effectively centered.
constants <- list(n.year=n.year,M=M,J=J,K1D=K1D,D.cov=data$D.cov,
                  n.cells=data$n.cells,n.cells.x=data$n.cells.x,
                  n.cells.y=data$n.cells.y,res=data$res,
                  x.vals=data$x.vals,y.vals=data$y.vals,
                  xlim=data$xlim,ylim=data$ylim,
                  cellArea=data$cellArea,n.cells=data$n.cells,
                  res=data$res)

#inits for Nimble
Niminits <- list(N=N.init,N.survive=N.survive.init,N.recruit=N.recruit.init,
                 ER=N.recruit.init,N.super=N.super.init,z.super=z.super.init,
                 z=z.init,z.start=z.start.init,z.stop=z.stop.init,
                 s=s.init,s.cell=s.cell.init,sigma.move=sigma.move.init,
                 beta0.phi=0,beta1.phi=0,D0=N.init[1]/(sum(InSS)*res^2),D.beta1=0,rsf.beta=0)

#data for Nimble
dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are
Nimdata <- list(y=y.nim,phi.cov=phi.cov.data,X=X.nim,dSS=dSS,
                dummy.data=dummy.data,cells=cells,InSS=data$InSS)

# set parameters to monitor
parameters <- c('N','gamma','N.recruit','N.survive','N.super','lambda.y1',
                'beta0.phi','beta1.phi','phi.cov.mu','phi.cov.sd','p0','sigma',
                'D0','D.beta1','rsf.beta','sigma.move')
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
               'phi.cov.mu','phi.cov.sd','p0','sigma','D0','D.beta1','rsf.beta','sigma.move')
# config.nodes <- c()
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

#activity center sampler. There are 3 samplers here for these cases
#1) z.super=1 and z=1, sSampler1 uses Metropolis-Hastings
#2) z.super=1 and z=0, sSampler2 uses Metropolis-Hastings with proposal sd tuned separately from above
#3) z.super=0, sSampler3 simulates entire new trajectory from priors.
#must identify and supply calcNodes in correct order.
for(i in 1:M){
  calcNodes2 <- c() #don't change anything below or this might not work
  for(g in 1:n.year){
    calcNodes1 <- Rmodel$getDependencies(paste0("s[",i,",",g,",1:2]"))
    if(g>1){ #must add node for s.cell in current year and s.lims for current year
      calcNodes.s.lims <- c(paste0("s.xlim[",i,",",g,",1:2]"),paste0("s.ylim[",i,",",g,",1:2]"))
      calcNodes1 <- c(paste0("s.cell[",i,", ",g,"]"),calcNodes.s.lims,calcNodes1)
    }
    conf$addSampler(target = paste0("s[",i,",",g,",1:2]"),
                    type = 'sSampler1',control=list(i=i,g=g,xlim=xlim,ylim=ylim,res=res,cells=cells,
                                                    calcNodes=calcNodes1,scale=1),silent = TRUE)
    conf$addSampler(target = paste0("s[",i,",",g,",1:2]"),
                    type = 'sSampler2',control=list(i=i,g=g,xlim=xlim,ylim=ylim,res=res,cells=cells,
                                                    calcNodes=calcNodes1,scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
  calcNodes2 <- c(paste0("s[",i,",1,1:2]"),paste0("s.cell[",i,",1]"),paste0("dummy.data[",i,"]")) #1st year nodes
  #subsequent year nodes
  for(g in 2:n.year){
    calcNodes2 <- c(calcNodes2,paste0("avail.dist[",i,",",g-1,",1:",n.cells,"]"),
                    paste0("use.dist[",i,",",g-1,",1:",n.cells,"]"),
                    paste0("s.cell[",i,",",g,"]"),
                    paste0("s.xlim[",i,",",g,",1:2]"),
                    paste0("s.ylim[",i,",",g,",1:2]"),
                    paste0("s[",i,",",g,",1:2]"))
  }
  conf$addSampler(target = paste0("s[",i,",1:",n.year,",1:2]"),
                  type = 'sSampler3',control=list(i=i,xlim=xlim,ylim=ylim,n.year=n.year,
                                                  n.cells=data$n.cells,n.cells.x=data$n.cells.x,
                                                  n.cells.y=data$n.cells.y,res=res,
                                                  x.vals=data$x.vals,y.vals=data$y.vals,
                                                  calcNodes=calcNodes2),silent = TRUE)
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
Cmcmc$run(5000,reset=FALSE) #can extend run by rerunning this line
end.time <- Sys.time()
time1 <- end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2 <- end.time-start.time2 # post-compilation run time

mvSamples <-  as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[-c(1:250),]))

#reminder what the targets are
data$N
data$N.recruit
data$N.survive
data$N[1]+sum(data$N.recruit) #N.super
