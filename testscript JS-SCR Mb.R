#This version has a trap-specific response to capture
#I am simulating and fitting year-specific p01 and p02 (first and subsequent baseline cap prob)
#you can change model file to include more structure
#I am fitting Mb the more efficient way by breaking
#the capture history into first and subsequent capture
#histories. So, no (secondary) occasion effects here.

library(nimble)
library(coda)
source("sim.JS.SCR.Mb.R")
source("Nimble Model JS-SCR Mb.R")
source("Nimble Functions JS-SCR Mb.R") #contains custom distributions and updates
source("sSampler Fixed.R") # activity center sampler that proposes from prior when z.super=0.
#this one works for fixed activity centers over years only

n.year <- 3 #number of years
lambda.y1 <- 100 #expected N in year 1
gamma <- rep(0.2,n.year-1) #yearly per-capita recruitment
beta0.phi <- qlogis(0.85) #survival intercept
beta1.phi <- 0.5 #phi response to individual covariate
p01 <- rep(0.15,n.year) #yearly detection probabilities at activity center, 1st capture (per trap)
p02 <- rep(0.5,n.year) #yearly detection probabilities at activity center, subsequent capture (per trap)
sigma <- rep(0.5,n.year) #yearly detection function scale
K <- rep(10,n.year) #yearly sampling occasions
buff <- 2 #state space buffer. Buffers maximal x and y dimensions of X below across years
X <- vector("list",n.year) #one trapping array per year
for(g in 1:n.year){ #using same trapping array every year here
  X[[g]] <- as.matrix(expand.grid(3:11,3:11))
}

data <- sim.JS.SCR.Mb(lambda.y1=lambda.y1,gamma=gamma,n.year=n.year,
            beta0.phi=beta0.phi,beta1.phi=beta1.phi,
            p01=p01, p02=p02,sigma=sigma,X=X,buff=buff,K=K)


##Initialize##
M <- 200 #data augmentation level. Check N.super posterior to make sure it never hits M
N.super.init <- nrow(data$y)
X <- data$X #pull X from data (won't be in environment if not simulated directly above)
if(N.super.init > M) stop("Must augment more than number of individuals captured")
J <- unlist(lapply(X,nrow)) #traps per year
J.max <- max(J)
K.max <- max(data$K)

#Mb stuff
#First! Let's split the n x n.year x J x K capture history into two of dimension
#n x n.year x J with one containing the number of first capture events and the second
#containing the number of subsequent capture events

y.tmp <- array(0,dim=c(M,n.year,J.max,K.max))
y.tmp[1:N.super.init,1:n.year,1:J.max,1:K.max] <- data$y #all these guys must be observed

#K3D (trap operation) is 1 if operational, 0 otherwise
#plugging in perfect operation here (use your own here if not perfect operation)
K3D <- array(NA,dim=c(n.year,J.max,K.max))
for(g in 1:n.year){
  K3D[g,1:J[g],1:K[g]] <- 1
}

n <- dim(data$y)[1]
y1 <- y2 <- array(0,dim=c(M,n.year,J.max))
K1 <- K2 <- array(0,dim=c(M,n.year,J.max))
for(g in 1:n.year){
  for(i in 1:M){
    for(j in 1:J[g]){
      position <- Position(y.tmp[i,g,j,],f=function(x){x>0})#occasion of first capture, NA if no capture
      if(is.na(position)){#no capture
        K1[i,g,j] <- sum(K3D[g,j,1:K[g]]) #sum trap op
        K2[i,g,j] <- 0 #no subsequent capture events
      }else{#there is a first capture
        K1[i,g,j] <- sum(K3D[g,j,1:position]) #sum trap op up to first capture event
        if(position<K[g]){#was first capture not last occasion?
          K2[i,g,j] <- sum(K3D[g,j,(position+1):K[g]]) #sum trap op after first capture event for subsequent capture
        }else{#otherwise, no subsequent capture events
          K2[i,g,j] <- 0
        }
        y1[i,g,j] <- y.tmp[i,g,j,position]
      }
    }
  }
}
y3D <- apply(y.tmp,c(1,2,3),sum)
y2 <- y3D-y1

#I have not tested that the algorithm above is correct in all cases of trap operation entered.
##Check for bugs here##
for(i in 1:M){
  for(g in 1:n.year){
    for(j in 1:J[g]){
      y.check <- y.tmp[i,g,j,]
      if(sum(y.check)>0){
        if(y1[i,g,j]!=1)stop("bug in y1")
        if(y2[i,g,j]!=(sum(y.check)-1))stop("bug in y2")
        first.cap <- which(y.check==1)[1]
        if(K1[i,g,j]!=first.cap)stop("bug in K1")
        if(K2[i,g,j]!=(sum(K3D[g,j,1:K[g]])-first.cap))stop("bug in K2")
      }else{
        if(y1[i,g,j]!=0)stop("bug in y1")
        if(y2[i,g,j]!=0)stop("bug in y2")
        if(K1[i,g,j]!=sum(K3D[g,j,1:K[g]]))stop("bug in K1")
        if(K2[i,g,j]!=0)stop("bug in K2")
      }
    }
  }
}

#does individual i have any exposures to subsequent captures in year y?
#currently not used, but could be used to determine which 
#individual-years subsequent capture likelihood needs to be calculated for
#could further specify by trap, but they're currently vector nodes. Not much more to be gained.
anyK2=apply(K2,c(1,2),function(x){any(x>0)})

#we are going to use y1, the first capture history with exposure K1 and y2, the subsequent
#capture history with exposure K2.

#initialize z, start with observed guys
y.nim <- apply(y.tmp,c(1,2,3),sum) #sum all captures over occasions
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
X.nim <- array(0,dim=c(n.year,J.max,2))
for(g in 1:n.year){
  X.nim[g,1:J[g],1:2] <- X[[g]]
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

#constants for Nimble
constants <- list(n.year=n.year, M=M, J=J, xlim=xlim, ylim=ylim, K1=K1, K2=K2)
#inits for Nimble
Niminits <- list(N=N.init,N.survive=N.survive.init,N.recruit=N.recruit.init,
                 ER=N.recruit.init,N.super=N.super.init,z.super=z.super.init,
                 s=s.init,
                 beta0.phi=beta0.phi,beta1.phi=beta1.phi,lambda.y1=lambda.y1) #for demonstration. Don't start at truth in practice.

#data for Nimble
Nimdata <- list(y1=y1,y2=y2,z=z.init,z.start=z.start.init,z.stop=z.stop.init,phi.cov=phi.cov.data,X=X.nim)

# set parameters to monitor
parameters <- c('N','gamma','N.recruit','N.survive','beta0.phi',
              'beta1.phi','lambda.y1','N.super','phi.cov.mu','phi.cov.sd',
              'p01','p02','sigma')
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
               'phi.cov.mu','phi.cov.sd','p01','p02','sigma')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,
                      nodes=config.nodes,useConjugacy = TRUE)

###*required* sampler replacements###
z.super.ups <- round(M*0.2) #how many z.super update proposals per iteration? 
#20% of M seems reasonable, but optimal will depend on data set
#loop here bc potentially different numbers of traps to vectorize in each year
y1.nodes <- y2.nodes <- pd1.nodes <- pd2.nodes <- kern.nodes <- c()
for(g in 1:n.year){
  y1.nodes <- c(y1.nodes,Rmodel$expandNodeNames(paste0("y1[1:",M,",",g,",1:",J[g],"]"))) #if you change y structure, change here
  y2.nodes <- c(y2.nodes,Rmodel$expandNodeNames(paste0("y2[1:",M,",",g,",1:",J[g],"]"))) #if you change y structure, change here
  pd1.nodes <- c(pd1.nodes,Rmodel$expandNodeNames(paste0("pd1[1:",M,",",g,",1:",J[g],"]"))) #if you change y structure, change here
  pd2.nodes <- c(pd2.nodes,Rmodel$expandNodeNames(paste0("pd2[1:",M,",",g,",1:",J[g],"]"))) #if you change y structure, change here
  kern.nodes <- c(kern.nodes,Rmodel$expandNodeNames(paste0("kern[1:",M,",",g,",1:",J[g],"]"))) #if you change y structure, change here
}
N.nodes <- Rmodel$expandNodeNames(paste0("N"))
N.survive.nodes <- Rmodel$expandNodeNames(paste0("N.survive[1:",n.year-1,"]"))
N.recruit.nodes <- Rmodel$expandNodeNames(paste0("N.recruit[1:",n.year-1,"]"))
ER.nodes <- Rmodel$expandNodeNames(paste0("ER[1:",n.year-1,"]"))
z.nodes <- Rmodel$expandNodeNames(paste0("z[1:",M,",1]"))
calcNodes <- c(N.nodes,N.recruit.nodes,y1.nodes,y2.nodes,z.nodes) #the ones that need likelihoods updated in mvSaved
conf$addSampler(target = c("z"),
                type = 'zSampler',control = list(M=M,n.year=n.year,z.obs=z.obs,z.super.ups=z.super.ups,
                                                 y1.nodes=y1.nodes,y2.nodes=y2.nodes,J=J,
                                                 pd1.nodes=pd1.nodes,pd2.nodes=pd2.nodes,
                                                 N.nodes=N.nodes,kern.nodes=kern.nodes,
                                                 z.nodes=z.nodes,ER.nodes=ER.nodes,
                                                 N.survive.nodes=N.survive.nodes,
                                                 N.recruit.nodes=N.recruit.nodes,
                                                 y2D=y.nim2D,calcNodes=calcNodes), silent = TRUE)

#activity center sampler. This sampler tunes activity centers when z.super[i]=1 and
#draws from the prior otherwise.
# conf$removeSampler(paste0("s[1:",M,", 1:2]")) #dont need to remove if not assigned by nimble
for(i in 1:M){
  conf$addSampler(target = paste0("s[",i,", 1:2]"),
                  type = 'sSampler',control=list(i=i,xlim=xlim,ylim=ylim,scale=1),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
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