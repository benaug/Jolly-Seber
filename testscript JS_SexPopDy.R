#This is a version with sex-specific population dynamics and recruitment
#I did not test this version, but this is the SCR version with space removed and I did
#test that one. No indications of any problems.

#Some notes:
#1) Dimensions probably prevent this from working with only 2 years of data (N.recruit, N.survive are not vectors)
#Need to modify custom updates in this case
#2) Object names that cannot be changed in the nimble model without changes in custom updates:
#N, N.recruit, N.survive, ER (male and female counterparts for all these, too),
# lambda.y1.M, lambda.y1.F, z.start, z.stop, z.obs
#phi[i,g] (must be of dimension M x n.year),
#Poisson assumptions on N.M[1], N.F[1] and N.recruit.M/N.recruit.F (but can include overdispersion with random effects)
#y can change dimension (e.g., for SCR), but need to account for that in defining "y.nodes"
#below to add custom updates. I think that is all...

library(nimble)
library(coda)
source("sim.JS.SexPopDy.R")
source("Nimble Model JS_SexPopDy.R")
source("Nimble Functions JS_SexPopDy.R") #contains custom distributions and updates


n.year <- 5 #number of years
lambda.y1.M <- 75 #expected male N in year 1
lambda.y1.F <- 125 #expected female N in year 1
#yearly per-capita recruitment
gamma.sex <- c(0.15,0.1) #male, then female, fixed across years
#sex-specific survival
phi.sex <- c(0.75,0.95) #male, then female, fixed across years
#sex-specific p
p.sex <- c(0.05,0.15) #male, then female, fixed across years
K <- rep(10,n.year) #yearly sampling occasions

#probability we observe sex for detected individuals (not a function of number of capture events)
#ASSUMPTION: sex observations are missing at random (same prob of observing male and female)
#model likely won't work well without most sexes observed
p.obs.sex <- 1

data <- sim.JS.SexPopDy(lambda.y1.M=lambda.y1.M,lambda.y1.F=lambda.y1.F,
                        n.year=n.year,gamma.sex=gamma.sex,phi.sex=phi.sex,
                        p.sex=p.sex,K=K,p.obs.sex=p.obs.sex)

##Initialize##
M <- 600 #data augmentation level.
# IMPORTANT: Check N.super posterior to make sure it never hits M. Otherwise, estimates will be biased.
N.super.init <- nrow(data$y)
if(N.super.init > M) stop("Must augment more than number of individuals captured")
y.nim <- matrix(0,M,n.year)
y.nim[1:N.super.init,] <- data$y #all these guys must be observed
#initialize z, start with observed guys
z.init <- 1*(y.nim>0)
z.init <- matrix(0,M,n.year)
z.start.init <- z.stop.init <- rep(NA,M)
for(i in 1:N.super.init){
  det.idx <- which(y.nim[i,]>0)
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

#now individual sex covariate. convert to 0 for female, 1 for male
sex.data <- c(data$sex,rep(NA,M-length(data$sex)))-1
sex.up <- which(is.na(sex.data)) #which individuals have missing cov values, used to determine which sexes need updating
#initialize missing sexes
sex.init <- sex.data
sex.init[sex.up] <- sample(c(0,1),length(sex.up),0.5)

#initialize N structures from z.init
N.M.init <- colSums(z.init[z.super.init==1&sex.init==0,])
N.F.init <- colSums(z.init[z.super.init==1&sex.init==1,])
N.init <- N.M.init + N.F.init
N.survive.M.init <- N.recruit.M.init <- rep(NA,n.year-1)
N.survive.F.init <- N.recruit.F.init <- rep(NA,n.year-1)
N.survive.init <- N.recruit.init <- rep(NA,n.year-1)
for(g in 2:n.year){
  N.survive.M.init[g-1] <- sum(z.init[,g-1]==1&z.init[,g]==1&z.super.init==1&sex.init==0)
  N.recruit.M.init[g-1] <- N.M.init[g]-N.survive.M.init[g-1]
  N.survive.F.init[g-1] <- sum(z.init[,g-1]==1&z.init[,g]==1&z.super.init==1&sex.init==1)
  N.recruit.F.init[g-1] <- N.F.init[g]-N.survive.F.init[g-1]
}
N.survive.init <- N.survive.M.init + N.survive.F.init
N.recruit.init <- N.recruit.M.init + N.recruit.F.init

#constants for Nimble
constants <- list(n.year=n.year, M=M, K=K)
#inits for Nimble
Niminits <- list(N=N.init,N.survive=N.survive.init,N.recruit=N.recruit.init,
                 N.M=N.M.init,N.survive.M=N.survive.M.init,N.recruit.M=N.recruit.M.init,
                 N.F=N.F.init,N.survive.F=N.survive.F.init,N.recruit.F=N.recruit.F.init,
                 ER.M=N.recruit.M.init,ER.F=N.recruit.F.init,
                 lambda.y1.M=N.M.init[1],lambda.y1.F=N.F.init[1],
                 N.super=N.super.init,z.super=z.super.init,
                 z=z.init,z.start=z.start.init,z.stop=z.stop.init,
                 sex=sex.init)

#data for Nimble
Nimdata <- list(y=y.nim,sex=sex.data)

# set parameters to monitor
parameters <- c('N','N.M',"N.F",'N.super',
                'N.recruit','N.recruit.M','N.recruit.F',
                'N.survive','N.survive.M','N.survive.F',
                'lambda.y1.M','lambda.y1.F','gamma.male','gamma.female','phi.sex',
                'p.sex')
parameters2 <- c('sex') #might want to monitor sex if interested in unobserved sex guy posteriors

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
config.nodes <- c('phi.sex','gamma.male','gamma.female','lambda.y1.M', #took sex out, put back in or make custom update
               'lambda.y1.F','p.sex')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,monitors2=parameters2,
                      nodes=config.nodes,useConjugacy = TRUE)

###*required* sampler replacements###
z.super.ups <- round(M*0.2) #how many z.super update proposals per iteration?
#20% of M seems reasonable, but optimal will depend on data set
y.nodes <- c(Rmodel$expandNodeNames(paste0("y[1:",M,",1:",n.year,"]")) )
N.nodes <- Rmodel$expandNodeNames(paste0("N"))
N.M.nodes <- Rmodel$expandNodeNames(paste0("N.M"))
N.F.nodes <- Rmodel$expandNodeNames(paste0("N.F"))
N.survive.nodes <- Rmodel$expandNodeNames(paste0("N.survive[1:",n.year-1,"]"))
N.survive.M.nodes <- Rmodel$expandNodeNames(paste0("N.survive.M[1:",n.year-1,"]"))
N.survive.F.nodes <- Rmodel$expandNodeNames(paste0("N.survive.F[1:",n.year-1,"]"))
N.recruit.nodes <- Rmodel$expandNodeNames(paste0("N.recruit[1:",n.year-1,"]"))
N.recruit.M.nodes <- Rmodel$expandNodeNames(paste0("N.recruit.M[1:",n.year-1,"]"))
N.recruit.F.nodes <- Rmodel$expandNodeNames(paste0("N.recruit.F[1:",n.year-1,"]"))
ER.M.nodes <- Rmodel$expandNodeNames(paste0("ER.M[1:",n.year-1,"]"))
ER.F.nodes <- Rmodel$expandNodeNames(paste0("ER.F[1:",n.year-1,"]"))
z.nodes <- Rmodel$expandNodeNames(paste0("z[1:",M,",1]"))
phi.nodes <-  Rmodel$expandNodeNames(paste0("phi"))

calcNodes <- c(N.nodes,N.recruit.nodes,
               N.M.nodes,N.recruit.M.nodes,
               N.F.nodes,N.recruit.F.nodes,
               N.survive.nodes,N.survive.M.nodes,N.survive.F.nodes,
               y.nodes,z.nodes,phi.nodes) #the ones that need likelihoods updated in mvSaved
conf$addSampler(target = c("z"),
                type = 'zSampler',control = list(M=M,n.year=n.year,
                                                 z.obs=z.obs,z.super.ups=z.super.ups,
                                                 y.nodes=y.nodes,
                                                 z.nodes=z.nodes,phi.nodes=phi.nodes,
                                                 ER.M.nodes=ER.M.nodes,
                                                 ER.F.nodes=ER.F.nodes,
                                                 N.nodes=N.nodes,N.M.nodes=N.M.nodes,N.F.nodes=N.F.nodes,
                                                 N.recruit.M.nodes=N.recruit.M.nodes,
                                                 N.recruit.F.nodes=N.recruit.F.nodes,
                                                 sex.up=sex.up,
                                                 calcNodes=calcNodes), silent = TRUE)


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
plot(mcmc(mvSamples[-c(1:50),]))

#reminder what the targets are
data$N
data$N.F
data$N.M
data$N.recruit
data$N.recruit.F
data$N.recruit.M
data$N.survive
data$N.survive.F
data$N.survive.M
data$N[1]+sum(data$N.recruit) #N.super


#if you record "sex" here
mvSamples2 <-  as.matrix(Cmcmc$mvSamples2)
latent.sex.idx <- sex.up[z.obs==1]
plot(mcmc(mvSamples2[-c(1:250),latent.sex.idx]))


# #Some sanity checks I used during debugging. Just checking that final
# #model states match between sex, z, and N structures
# 
# #check N
# N.count <- rep(NA,n.year)
# N.count.M <- rep(NA,n.year)
# N.count.F <- rep(NA,n.year)
# for(g2 in 1:n.year){
#   N.count[g2] <- sum(Cmodel$z[Cmodel$z.super==1,g2]==1)
#   N.count.M[g2] <- sum(Cmodel$z[Cmodel$z.super==1,g2]==1&Cmodel$sex[Cmodel$z.super==1]==0)
#   N.count.F[g2] <- sum(Cmodel$z[Cmodel$z.super==1,g2]==1&Cmodel$sex[Cmodel$z.super==1]==1)
# }
# all(N.count==Cmodel$N)
# all(N.count.M==Cmodel$N.M)
# all(N.count.F==Cmodel$N.F)
# 
# 
# #check N.recruit count
# N.count <- rep(NA,n.year-1)
# N.count.M <- rep(NA,n.year-1)
# N.count.F <- rep(NA,n.year-1)
# for(g2 in 2:n.year){
#   N.count[g2-1] <- sum(Cmodel$z[Cmodel$z.super==1,g2-1]==0&Cmodel$z[Cmodel$z.super==1,g2]==1)
#   N.count.M[g2-1] <- sum(Cmodel$z[Cmodel$z.super==1,g2-1]==0&Cmodel$z[Cmodel$z.super==1,g2]==1&Cmodel$sex[Cmodel$z.super==1]==0)
#   N.count.F[g2-1] <- sum(Cmodel$z[Cmodel$z.super==1,g2-1]==0&Cmodel$z[Cmodel$z.super==1,g2]==1&Cmodel$sex[Cmodel$z.super==1]==1)
# }
# all(N.count==Cmodel$N.recruit)
# all(N.count.M==Cmodel$N.recruit.M)
# all(N.count.F==Cmodel$N.recruit.F)
# 
# #check N.survive count
# N.count <- rep(NA,n.year-1)
# N.count.M <- rep(NA,n.year-1)
# N.count.F <- rep(NA,n.year-1)
# for(g2 in 2:n.year){
#   N.count[g2-1] <- sum(Cmodel$z[Cmodel$z.super==1,g2-1]==1&Cmodel$z[Cmodel$z.super==1,g2]==1)
#   N.count.M[g2-1] <- sum(Cmodel$z[Cmodel$z.super==1,g2-1]==1&Cmodel$z[Cmodel$z.super==1,g2]==1&Cmodel$sex[Cmodel$z.super==1]==0)
#   N.count.F[g2-1] <- sum(Cmodel$z[Cmodel$z.super==1,g2-1]==1&Cmodel$z[Cmodel$z.super==1,g2]==1&Cmodel$sex[Cmodel$z.super==1]==1)
# }
# all(N.count==Cmodel$N.survive)
# all(N.count.M==Cmodel$N.survive.M)
# all(N.count.F==Cmodel$N.survive.F)
# 
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

