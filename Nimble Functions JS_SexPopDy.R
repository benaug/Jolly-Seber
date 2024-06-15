#this is used to restrict likelihood evaluation to only the years relevant for survival for each individual
dSurvival <- nimbleFunction(
  run = function(x = double(1), phi = double(1), z.start = double(0), z.stop = double(0),
                 log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    n.year <- length(phi)+1
    #extract first and last survival event years
    surv.start <- z.start+1
    surv.stop <- z.stop+1 #count death events, first z[i,]=0
    if(surv.start <= n.year){ #if surv.start beyond last year, no survival events, logProb=0
      if(surv.stop > n.year){ #but can't survive past n.year
        surv.stop <- n.year 
      }
      for(g in surv.start:surv.stop){ #sum logprob over survival event years
        logProb <- logProb + dbinom(x[g], size = 1, p = phi[g-1], log = TRUE)
      }
    }
    return(logProb)
  }
)

#make dummy random vector generator to make nimble happy
rSurvival <- nimbleFunction(
  run = function(n = integer(0),phi = double(1), z.start = double(0), z.stop = double(0)) {
    returnType(double(1))
    n.year <- length(phi)
    return(rep(0,n.year))
  }
)

#custom observation model distribution that makes custom updates easier
dbinomial2 <- nimbleFunction(
  run = function(x = double(0), p = double(0), K = double(0), z = double(0), z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z.super==0 | z.super==1&z==0){#skip calculation if not is superpop, or in superpop, but not alive in this year
      return(0)
    }else{
      logProb <- dbinom(x, size = K, p = p, log = TRUE)
      return(logProb)
    }
  }
)

rbinomial2 <- nimbleFunction(
  run = function(n = integer(0), p = double(0), K = double(0), z = double(0), z.super = double(0)) {
    returnType(double(0))
    return(0)
  }
)


#all z updates live here
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    z.super.ups <- control$z.super.ups
    n.year <- control$n.year
    z.obs <- control$z.obs
    z.nodes <- control$z.nodes
    y.nodes <- control$y.nodes
    phi.nodes <- control$phi.nodes
    N.nodes <- control$N.nodes
    N.M.nodes <- control$N.M.nodes
    N.F.nodes <- control$N.F.nodes
    ER.M.nodes <- control$ER.M.nodes
    ER.F.nodes <- control$ER.F.nodes
    N.recruit.M.nodes <- control$N.recruit.M.nodes
    N.recruit.F.nodes <- control$N.recruit.F.nodes
    sex.up <- control$sex.up
    calcNodes <- control$calcNodes
  },
  run = function() {
    # 1) Detected guy updates: z.start, z.stop
    # 1a) z start update (z.stop update below)
    for(i in 1:M){
      if(z.obs[i]==1&model$y[i,1]==0){ #for detected guys, skip if observed 1st year
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        N.curr <- model$N
        N.recruit.curr <- model$N.recruit
        if(model$sex[i]==0){
          N.M.curr <- model$N.M
          N.recruit.M.curr <- model$N.recruit.M
        }else{
          N.F.curr <- model$N.F
          N.recruit.F.curr <- model$N.recruit.F
        }
        y <- model$y[i,]
        dets <- which(model$y[i,]>0)
        first.det <- min(dets)
        lp.start <- rep(-Inf,n.year)
        i.idx <- seq(i,M*n.year,M) #used to reference correct y nodes
        for(g in 1:first.det){ #must be recruited in year with first detection or before
          z.start.prop <- g
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.year)
          z.prop[g:first.det] <- 1 #must be alive until first detection
          if(first.det < n.year){
            z.prop[(first.det+1):n.year] <- z.curr[(first.det+1):n.year] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop

          #update N, N.recruit, N.survive. These individuals always in superpopulation
          #1) Update N
          model$N <<- N.curr - z.curr + z.prop
          #2) Update N.recruit
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- N.recruit.curr[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g

          #now repeat for sex
          if(model$sex[i]==0){ #male
            #1) Update N
            model$N.M <<- N.M.curr - z.curr + z.prop
            #2) Update N.recruit
            model$N.recruit.M <<- N.recruit.M.curr #set back to original first
            if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
              model$N.recruit.M[z.start.curr-1] <<- N.recruit.M.curr[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
              model$N.recruit.M[z.start.prop-1] <<- N.recruit.M.curr[z.start.prop-1] + 1
            }
            #3) Update N.survive
            model$N.survive.M <<- model$N.M[2:n.year]-model$N.recruit.M #survivors are guys alive in year g-1 minus recruits in this year g
          }else{ #female
            #1) Update N
            model$N.F <<- N.F.curr - z.curr + z.prop
            #2) Update N.recruit
            model$N.recruit.F <<- N.recruit.F.curr #set back to original first
            if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
              model$N.recruit.F[z.start.curr-1] <<- N.recruit.F.curr[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
              model$N.recruit.F[z.start.prop-1] <<- N.recruit.F.curr[z.start.prop-1] + 1
            }
            #3) Update N.survive
            model$N.survive.F <<- model$N.F[2:n.year]-model$N.recruit.F #survivors are guys alive in year g-1 minus recruits in this year g
          }
          # recruit likelihood conditional on having recruited (or alive in year 1)
          #must account for sex-specificity, Updating N changes both ER.M and ER.F
          model$calculate(ER.M.nodes)
          model$calculate(ER.F.nodes)
          recruit.probs.M <- c(model$lambda.y1.M,model$ER.M)
          recruit.probs.F <- c(model$lambda.y1.F,model$ER.F)
          recruit.probs.M <- recruit.probs.M/sum(recruit.probs.M)
          recruit.probs.F <- recruit.probs.F/sum(recruit.probs.F)
          lp.recruit <- 0
          for(i2 in 1:M){
            if(model$sex[i2]==0){
              lp.recruit <- lp.recruit + log(recruit.probs.M[model$z.start[i2]])
            }else{
              lp.recruit <- lp.recruit + log(recruit.probs.F[model$z.start[i2]])
            }
          }
          lp.y <- model$calculate(y.nodes[i.idx])
          lp.surv <- model$calculate(z.nodes[i])
          lp.start[g] <- lp.recruit + lp.y + lp.surv
        }
        maxlp <- max(lp.start) #deal with overflow
        prop.probs <- exp(lp.start-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)

        z.start.prop <- rcat(1,prop.probs)
        model$z.start[i] <<- z.start.curr #set back to original

        if(model$z.start[i]!=z.start.prop){#if proposal is same as current, no need to replace anything
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.year)
          z.prop[model$z.start[i]:first.det] <- 1 #must be alive until first detection
          if(first.det < n.year){
            z.prop[(first.det+1):n.year] <- z.curr[(first.det+1):n.year] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- N.recruit.curr[z.start.prop-1] + 1
          }
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          #now repeat for sex
          if(model$sex[i]==0){ #male
            model$N.M <<- N.M.curr - z.curr + z.prop
            model$N.recruit.M <<- N.recruit.M.curr #set back to original first
            if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
              model$N.recruit.M[z.start.curr-1] <<- N.recruit.M.curr[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
              model$N.recruit.M[z.start.prop-1] <<- N.recruit.M.curr[z.start.prop-1] + 1
            }
            model$N.survive.M <<- model$N.M[2:n.year]-model$N.recruit.M #survivors are guys alive in year g-1 minus recruits in this year g
          }else{ #female
            model$N.F <<- N.F.curr - z.curr + z.prop
            model$N.recruit.F <<- N.recruit.F.curr #set back to original first
            if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
              model$N.recruit.F[z.start.curr-1] <<- N.recruit.F.curr[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
              model$N.recruit.F[z.start.prop-1] <<- N.recruit.F.curr[z.start.prop-1] + 1
            }
            model$N.survive.F <<- model$N.F[2:n.year]-model$N.recruit.F #survivors are guys alive in year g-1 minus recruits in this year g
          }
          model$calculate(ER.M.nodes)
          model$calculate(ER.F.nodes)
          #update these logProbs
          model$calculate(y.nodes[i.idx])
          model$calculate(z.nodes[i])
          model$calculate(N.M.nodes[1])
          model$calculate(N.F.nodes[1])
          model$calculate(N.recruit.M.nodes)
          model$calculate(N.recruit.F.nodes)
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          if(model$sex[i]==0){ #male
            mvSaved["N.M",1] <<- model[["N.M"]]
            mvSaved["N.survive.M",1] <<- model[["N.survive.M"]]
            mvSaved["N.recruit.M",1] <<- model[["N.recruit.M"]]
          }else{ #female
            mvSaved["N.F",1] <<- model[["N.F"]]
            mvSaved["N.survive.F",1] <<- model[["N.survive.F"]]
            mvSaved["N.recruit.F",1] <<- model[["N.recruit.F"]]
          }
          mvSaved["ER.M",1] <<- model[["ER.M"]]
          mvSaved["ER.F",1] <<- model[["ER.F"]]
        }else{
          model[["z.start"]][i] <<- mvSaved["z.start",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["N.recruit"]] <<- mvSaved["N.recruit",1]
          if(model$sex[i]==0){ #male
            model[["N.M"]] <<- mvSaved["N.M",1]
            model[["N.survive.M"]] <<- mvSaved["N.survive.M",1]
            model[["N.recruit.M"]] <<- mvSaved["N.recruit.M",1]
          }else{ #female
            model[["N.F"]] <<- mvSaved["N.F",1]
            model[["N.survive.F"]] <<- mvSaved["N.survive.F",1]
            model[["N.recruit.F"]] <<- mvSaved["N.recruit.F",1]
          }
          model[["ER.M"]] <<- mvSaved["ER.M",1]
          model[["ER.F"]] <<- mvSaved["ER.F",1]
          #set these logProbs back
          model$calculate(y.nodes[i.idx])
          model$calculate(N.M.nodes[1])
          model$calculate(N.F.nodes[1])
          model$calculate(N.recruit.M.nodes)
          model$calculate(N.recruit.F.nodes)
          model$calculate(z.nodes[i])
        }
      }
    }

    #1b) z stop update (z.start update above)
    for(i in 1:M){
      if(z.obs[i]==1&model$y[i,n.year]==0){ #for detected guys, skip if observed in final year
        z.curr <- model$z[i,]
        z.stop.curr <- model$z.stop[i]
        N.curr <- model$N
        if(model$sex[i]==0){
          N.M.curr <- model$N.M
        }else{
          N.F.curr <- model$N.F
        }
        y <- model$y[i,]
        dets <- which(model$y[i,]>0)
        last.det <- max(dets)
        lp.stop <- rep(-Inf,n.year)
        i.idx <- seq(i,M*n.year,M) #used to reference correct y nodes
        for(g in (last.det):n.year){ #can't die on or before year of last detection
          model$z.stop[i] <<- g
          z.prop <- rep(0,n.year)
          z.prop[last.det:g] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          #update N, number of recruits does not change going backwards
          model$N <<- N.curr - z.curr + z.prop
          if(model$sex[i]==0){ #male
            model$N.M <<- N.M.curr - z.curr + z.prop
          }else{ #female
            model$N.F <<- N.F.curr - z.curr + z.prop
          }
          # recruit likelihood conditional on having recruited (or alive in year 1)
          #must account for sex-specificity, Updating N changes both ER.M and ER.F
          model$calculate(ER.M.nodes)
          model$calculate(ER.F.nodes)
          recruit.probs.M <- c(model$lambda.y1.M,model$ER.M)
          recruit.probs.F <- c(model$lambda.y1.F,model$ER.F)
          recruit.probs.M <- recruit.probs.M/sum(recruit.probs.M)
          recruit.probs.F <- recruit.probs.F/sum(recruit.probs.F)
          lp.recruit <- 0
          for(i2 in 1:M){
            if(model$sex[i2]==0){
              lp.recruit <- lp.recruit + log(recruit.probs.M[model$z.start[i2]])
            }else{
              lp.recruit <- lp.recruit + log(recruit.probs.F[model$z.start[i2]])
            }
          }
          lp.y <- model$calculate(y.nodes[i.idx])
          lp.surv <- model$calculate(z.nodes[i])
          lp.stop[g] <- lp.recruit + lp.y + lp.surv
        }
        maxlp <- max(lp.stop) #deal with overflow
        prop.probs <- exp(lp.stop-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        z.stop.prop <- rcat(1,prop.probs)
        model$z.stop[i] <<- z.stop.curr #set back to original
        if(model$z.stop[i]!=z.stop.prop){#if proposal differs from current
          model$z.stop[i] <<- z.stop.prop
          z.prop <- rep(0,n.year)
          z.prop[last.det:model$z.stop[i]] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          if(model$sex[i]==0){ #male
            model$N.M <<- N.M.curr - z.curr + z.prop
            model$N.survive.M <<- model$N.M[2:n.year]-model$N.recruit.M
          }else{ #female
            model$N.F <<- N.F.curr - z.curr + z.prop
            model$N.survive.F <<- model$N.F[2:n.year]-model$N.recruit.F
          }
          model$calculate(ER.M.nodes)
          model$calculate(ER.F.nodes)
          #update these logProbs
          model$calculate(N.M.nodes[1])
          model$calculate(N.F.nodes[1])
          model$calculate(N.recruit.M.nodes)
          model$calculate(N.recruit.F.nodes)
          model$calculate(y.nodes[i.idx])
          model$calculate(z.nodes[i])
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          if(model$sex[i]==0){ #male
            mvSaved["N.M",1] <<- model[["N.M"]]
            mvSaved["N.survive.M",1] <<- model[["N.survive.M"]]
          }else{ #female
            mvSaved["N.F",1] <<- model[["N.F"]]
            mvSaved["N.survive.F",1] <<- model[["N.survive.F"]]
          }
          mvSaved["ER.M",1] <<- model[["ER.M"]]
          mvSaved["ER.F",1] <<- model[["ER.F"]]
        }else{
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          if(model$sex[i]==0){ #male
            model[["N.M"]] <<- mvSaved["N.M",1]
            model[["N.survive.M"]] <<- mvSaved["N.survive.M",1]
          }else{ #female
            model[["N.F"]] <<- mvSaved["N.F",1]
            model[["N.survive.F"]] <<- mvSaved["N.survive.F",1]
          }
          model[["ER.M"]] <<- mvSaved["ER.M",1]
          model[["ER.F"]] <<- mvSaved["ER.F",1]
          #set these logProbs back
          model$calculate(y.nodes[i.idx])
          model$calculate(N.M.nodes[1])
          model$calculate(N.F.nodes[1])
          model$calculate(N.recruit.M.nodes)
          model$calculate(N.recruit.F.nodes)
          model$calculate(z.nodes[i])
        }
      }
    }
    #2) undetected guy update. Regardless of whether they are in or out of superpopulation.
    #simulate new sex + recruitment/survival
    for(i in 1:M){
      if(z.obs[i]==0){
        initial.sex=model$sex[i] #store this for use below
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        z.stop.curr <- model$z.stop[i]
        i.idx <- seq(i,M*n.year,M) #used to reference correct y nodes
        i.idx2 <- seq(i,M*(n.year-1),M) #used to reference correct phi nodes

        #paste male and female recruit probs, Poisson RVs conditioned on total (N.super)
        recruit.probs.curr <- c(model$lambda.y1.M,model$ER.M,model$lambda.y1.F,model$ER.F)
        recruit.probs.curr <- recruit.probs.curr/sum(recruit.probs.curr)
        recruit.idx.curr <- c(model$z.start+model$sex*n.year) #used to reference pasted recruit.probs

        #record backwards phi proposal probs
        log.prop.back.surv <- 0
        if(z.start.curr < n.year){#if you don't recruit in final year
          for(g in (z.start.curr+1):n.year){
            log.prop.back.surv <- log.prop.back.surv + dbinom(z.curr[g],1,model$phi[i,g-1]*z.curr[g-1],log=TRUE)
          }
        }

        #Get Initial LogProbs
        # lp.initial.sex <- model$getLogProb(sex.nodes[i])
        lp.initial.surv <- model$getLogProb(z.nodes[i])
        lp.initial.y <- model$getLogProb(y.nodes[i.idx])
        if(model$z.super[i]==1){ #need to consider all individuals
          lp.initial.recruit <- sum(log(recruit.probs.curr[recruit.idx.curr]))
        }else{ #only need to consider focal individual
          lp.initial.recruit = log(recruit.probs.curr[recruit.idx.curr[i]])
        }

        #start updates#

        #simulate recruitment, update z.start and sex
        recruit.prop <- rcat(1,recruit.probs.curr)
        z.prop <- rep(0,n.year)
        log.prop.for.recruit <- log(recruit.probs.curr[recruit.prop])
        if(recruit.prop<=n.year){ #simulated male
          z.start.prop <- recruit.prop
          model$sex[i] <<- 0
        }else{ #simulated female
          z.start.prop <- recruit.prop - n.year
          model$sex[i] <<- 1
        }
        z.prop[z.start.prop] <- 1
        #update phi bc sex can change
        model$calculate(phi.nodes[i.idx2])
        #recruit.idx can change
        recruit.idx.prop <- recruit.idx.curr
        recruit.idx.prop[i] <- recruit.prop

        #simulate survival with updated phi
        log.prop.for.surv <- 0
        if(z.start.prop < n.year){#if you don't recruit in final year
          for(g in (z.start.prop+1):n.year){
            z.prop[g] <- rbinom(1,1,model$phi[i,g-1]*z.prop[g-1])
            log.prop.for.surv <- log.prop.for.surv + dbinom(z.prop[g],1,model$phi[i,g-1]*z.prop[g-1],log=TRUE)
          }
        }
        z.on.prop <- which(z.prop==1)
        z.stop.prop <- max(z.on.prop)
        model$z[i,] <<- z.prop
        model$z.start[i] <<- z.start.prop
        model$z.stop[i] <<- z.stop.prop

        #update N, N.recruit, N.survive only if individual is in superpopulation
        if(model$z.super[i]==1){
          #1) Update N
          model$N <<- model$N - z.curr + z.prop
          #2) Update N.recruit
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          #repeat for sex
          if(initial.sex==0&model$sex[i]==0){ #male to male
            model$N.M <<- model$N.M - z.curr + z.prop
            if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
              model$N.recruit.M[z.start.curr-1] <<- model$N.recruit.M[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
              model$N.recruit.M[z.start.prop-1] <<- model$N.recruit.M[z.start.prop-1] + 1
            }
            model$N.survive.M <<- model$N.M[2:n.year]-model$N.recruit.M #survivors are guys alive in year g-1 minus recruits in this year g
          }else if(initial.sex==1&model$sex[i]==1){ #female to female
            model$N.F <<- model$N.F - z.curr + z.prop
            if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
              model$N.recruit.F[z.start.curr-1] <<- model$N.recruit.F[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
              model$N.recruit.F[z.start.prop-1] <<- model$N.recruit.F[z.start.prop-1] + 1
            }
            model$N.survive.F <<- model$N.F[2:n.year]-model$N.recruit.F #survivors are guys alive in year g-1 minus recruits in this year g
          }else if(initial.sex==0&model$sex[i]==1){ #male to female
            #subtract current z from males, add new z to females
            model$N.M <<- model$N.M - z.curr
            model$N.F <<- model$N.F + z.prop
            if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
              model$N.recruit.M[z.start.curr-1] <<- model$N.recruit.M[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
              model$N.recruit.F[z.start.prop-1] <<- model$N.recruit.F[z.start.prop-1] + 1
            }
            model$N.survive.M <<- model$N.M[2:n.year]-model$N.recruit.M #survivors are guys alive in year g-1 minus recruits in this year g
            model$N.survive.F <<- model$N.F[2:n.year]-model$N.recruit.F #survivors are guys alive in year g-1 minus recruits in this year g

          }else if(initial.sex==1&model$sex[i]==0){ #female to male
            #subtract current z from females, add new z to males
            model$N.F <<- model$N.F - z.curr
            model$N.M <<- model$N.M + z.prop
            if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
              model$N.recruit.F[z.start.curr-1] <<- model$N.recruit.F[z.start.curr-1] - 1
            }
            if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
              model$N.recruit.M[z.start.prop-1] <<- model$N.recruit.M[z.start.prop-1] + 1
            }
            model$N.survive.M <<- model$N.M[2:n.year]-model$N.recruit.M #survivors are guys alive in year g-1 minus recruits in this year g
            model$N.survive.F <<- model$N.F[2:n.year]-model$N.recruit.F #survivors are guys alive in year g-1 minus recruits in this year g
          }
          model$calculate(ER.M.nodes) #update ER when N updated
          model$calculate(ER.F.nodes) #update ER when N updated
          # Compute new recruit.probs when N updated
          recruit.probs.prop <- c(model$lambda.y1.M,model$ER.M,model$lambda.y1.F,model$ER.F)
          recruit.probs.prop <- recruit.probs.prop/sum(recruit.probs.prop)
        }else{
          recruit.probs.prop <- recruit.probs.curr
        }

        #record recruit backprobs
        log.prop.back.recruit <- log(recruit.probs.prop[recruit.idx.curr[i]])

        if(model$z.super[i]==1){ #need to consider all individuals
          lp.proposed.recruit <- sum(log(recruit.probs.prop[recruit.idx.prop]))
        }else{ #only need to consider focal individual
          lp.proposed.recruit = log(recruit.probs.prop[recruit.idx.prop[i]])
        }

        #update pd nodes when z/sex changes
        lp.proposed.y <- model$calculate(y.nodes[i.idx])
        lp.proposed.surv <- model$calculate(z.nodes[i])

        #Add up likelihoods and prop probs
        #sex and survival likelihoods and prop probs cancel out, but leaving them in here
        #recruitment likelihood and prop probs cancel when z.super[i]=0, leaving in
        lp.initial.total <- lp.initial.recruit + lp.initial.surv + lp.initial.y
        lp.proposed.total <- lp.proposed.recruit + lp.proposed.surv + lp.proposed.y

        if(model$z.super[i]==1){ #survival likelihoods and prop probs cancel, nothing else...
          log.prop.for <- log.prop.for.recruit + log.prop.for.surv
          log.prop.back <-log.prop.back.recruit + log.prop.back.surv
        }else{ #log_MH_ratio always 0 if z.super[i]=0. Calculating lots of extra stuff not required...
          log.prop.for <-  log.prop.for.recruit + log.prop.for.surv
          log.prop.back <- log.prop.back.recruit + log.prop.back.surv
        }

        #MH step
        log_MH_ratio <- (lp.proposed.total + log.prop.back) - (lp.initial.total + log.prop.for)
        accept <- decide(log_MH_ratio)
        if(accept){
          model$calculate(N.M.nodes[1]) #update if accepted
          model$calculate(N.F.nodes[1]) #update if accepted
          model$calculate(N.recruit.M.nodes)
          model$calculate(N.recruit.F.nodes)
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["sex",1][i] <<- model[["sex"]][i]
          if(model$z.super[i]==1){
            mvSaved["ER.M",1] <<- model[["ER.M"]]
            mvSaved["ER.F",1] <<- model[["ER.F"]]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.M",1] <<- model[["N.M"]]
            mvSaved["N.survive.M",1] <<- model[["N.survive.M"]]
            mvSaved["N.recruit.M",1] <<- model[["N.recruit.M"]]
            mvSaved["N.F",1] <<- model[["N.F"]]
            mvSaved["N.survive.F",1] <<- model[["N.survive.F"]]
            mvSaved["N.recruit.F",1] <<- model[["N.recruit.F"]]
          }
        }else{
          model[["z.start"]][i] <<- mvSaved["z.start",1][i]
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["sex"]][i] <<- mvSaved["sex",1][i]
          if(model$z.super[i]==1){
            model[["ER.M"]] <<- mvSaved["ER.M",1]
            model[["ER.F"]] <<- mvSaved["ER.F",1]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.M"]] <<- mvSaved["N.M",1]
            model[["N.survive.M"]] <<- mvSaved["N.survive.M",1]
            model[["N.recruit.M"]] <<- mvSaved["N.recruit.M",1]
            model[["N.F"]] <<- mvSaved["N.F",1]
            model[["N.survive.F"]] <<- mvSaved["N.survive.F",1]
            model[["N.recruit.F"]] <<- mvSaved["N.recruit.F",1]
          }
          #set these logProbs back
          model$calculate(N.M.nodes[1])
          model$calculate(N.F.nodes[1])
          model$calculate(N.recruit.M.nodes)
          model$calculate(N.recruit.F.nodes)
          model$calculate(phi.nodes[i.idx2])
          model$calculate(z.nodes[i])
          model$calculate(y.nodes[i.idx])
        }
      }
    }

    #3) update z.super
    for(up in 1:z.super.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected individual
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z.super==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        if(z.obs[pick]==1){ #is this individual detected?
          reject <- TRUE #if so, we reject (could never select these inds, but then need to account for asymmetric proposal)
        }
        if(!reject){
          pick.idx <- seq(pick,M*n.year,M) #used to reference correct y nodes
          #get initial logprobs
          lp.initial.N.M <- model$getLogProb(N.M.nodes[1])
          lp.initial.N.F <- model$getLogProb(N.F.nodes[1])
          lp.initial.N.recruit.M <- model$getLogProb(N.recruit.M.nodes)
          lp.initial.N.recruit.F <- model$getLogProb(N.recruit.F.nodes)
          lp.initial.y <- model$getLogProb(y.nodes[pick.idx])

          # propose new N.super/z.super
          model$N.super <<-  model$N.super - 1
          model$z.super[pick] <<- 0

          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N - model$z[pick,]
          #2) Update N.recruit
          if(model$z.start[pick] > 1){ #if wasn't in pop in year 1
            model$N.recruit[model$z.start[pick]-1] <<- model$N.recruit[model$z.start[pick]-1] - 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.year]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          #repeat for sex
          if(model$sex[pick]==0){
            model$N.M <<- model$N.M - model$z[pick,]
            if(model$z.start[pick] > 1){ #if wasn't in pop in year 1
              model$N.recruit.M[model$z.start[pick]-1] <<- model$N.recruit.M[model$z.start[pick]-1] - 1
            }
            #3) Update N.survive
            model$N.survive.M <<- model$N.M[2:n.year]-model$N.recruit.M #survivors are guys alive in year g-1 minus recruits in this year g
          }else{
            model$N.F <<- model$N.F - model$z[pick,]
            if(model$z.start[pick] > 1){ #if wasn't in pop in year 1
              model$N.recruit.F[model$z.start[pick]-1] <<- model$N.recruit.F[model$z.start[pick]-1] - 1
            }
            #3) Update N.survive
            model$N.survive.F <<- model$N.F[2:n.year]-model$N.recruit.F #survivors are guys alive in year g-1 minus recruits in this year g
          }

          model$calculate(ER.M.nodes) #update ER when N updated
          model$calculate(ER.F.nodes) #update ER when N updated
          #get proposed logprobs for N and y
          lp.proposed.N.M <- model$calculate(N.M.nodes[1])
          lp.proposed.N.F <- model$calculate(N.F.nodes[1])
          lp.proposed.N.recruit.M <- model$calculate(N.recruit.M.nodes)
          lp.proposed.N.recruit.F <- model$calculate(N.recruit.F.nodes)
          lp.proposed.y <- model$calculate(y.nodes[pick.idx]) #will always be 0

          lp.initial.total <- lp.initial.N.M + lp.initial.N.F + lp.initial.y + lp.initial.N.recruit.M + lp.initial.N.recruit.F
          lp.proposed.total <- lp.proposed.N.M + lp.proposed.N.F + lp.proposed.y + lp.proposed.N.recruit.M + lp.proposed.N.recruit.F

          #MH step
          log_MH_ratio <- lp.proposed.total - lp.initial.total
          accept <- decide(log_MH_ratio)

          if(accept) {
            mvSaved["z.super",1] <<- model[["z.super"]]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            if(model$sex[pick]==0){
              mvSaved["N.M",1] <<- model[["N.M"]]
              mvSaved["N.survive.M",1] <<- model[["N.survive.M"]]
              mvSaved["N.recruit.M",1] <<- model[["N.recruit.M"]]
            }else{
              mvSaved["N.F",1] <<- model[["N.F"]]
              mvSaved["N.survive.F",1] <<- model[["N.survive.F"]]
              mvSaved["N.recruit.F",1] <<- model[["N.recruit.F"]]
            }
            mvSaved["ER.M",1] <<- model[["ER.M"]]
            mvSaved["ER.F",1] <<- model[["ER.F"]]
          }else{
            model[["z.super"]] <<- mvSaved["z.super",1]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            if(model$sex[pick]==0){
              model[["N.M"]] <<- mvSaved["N.M",1]
              model[["N.survive.M"]] <<- mvSaved["N.survive.M",1]
              model[["N.recruit.M"]] <<- mvSaved["N.recruit.M",1]
            }else{
              model[["N.F"]] <<- mvSaved["N.F",1]
              model[["N.survive.F"]] <<- mvSaved["N.survive.F",1]
              model[["N.recruit.F"]] <<- mvSaved["N.recruit.F",1]
            }
            model[["ER.M"]] <<- mvSaved["ER.M",1]
            model[["ER.F"]] <<- mvSaved["ER.F",1]
            #set these logProbs back
            model$calculate(N.M.nodes[1])
            model$calculate(N.F.nodes[1])
            model$calculate(y.nodes[pick.idx])
            model$calculate(N.recruit.M.nodes)
            model$calculate(N.recruit.F.nodes)
          }
        }
      }else{#add
        if(model$N.super[1] < M){ #cannot update if z.super maxed out. Need to raise M
          z.off <- which(model$z.super==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          pick.idx <- seq(pick,M*n.year,M)

          #get initial logProbs
          lp.initial.N.M <- model$getLogProb(N.M.nodes[1])
          lp.initial.N.F <- model$getLogProb(N.F.nodes[1])
          lp.initial.N.recruit.M <- model$getLogProb(N.recruit.M.nodes)
          lp.initial.N.recruit.F <- model$getLogProb(N.recruit.F.nodes)
          lp.initial.y <- model$getLogProb(y.nodes[pick.idx]) #will always be 0

          #propose new N/z
          model$N.super <<-  model$N.super + 1
          model$z.super[pick] <<- 1

          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N + model$z[pick,]
          #2) Update N.recruit
          if(model$z.start[pick] > 1){ #if wasn't in pop in year 1
            model$N.recruit[model$z.start[pick]-1] <<- model$N.recruit[model$z.start[pick]-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.year] - model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          #repeat for sex
          if(model$sex[pick]==0){
            model$N.M <<- model$N.M + model$z[pick,]
            if(model$z.start[pick] > 1){ #if wasn't in pop in year 1
              model$N.recruit.M[model$z.start[pick]-1] <<- model$N.recruit.M[model$z.start[pick]-1] + 1
            }
            #3) Update N.survive
            model$N.survive.M <<- model$N.M[2:n.year] - model$N.recruit.M #survivors are guys alive in year g-1 minus recruits in this year g
            model$calculate(ER.M.nodes) #update ER when N updated
          }else{
            model$N.F <<- model$N.F + model$z[pick,]
            if(model$z.start[pick] > 1){ #if wasn't in pop in year 1
              model$N.recruit.F[model$z.start[pick]-1] <<- model$N.recruit.F[model$z.start[pick]-1] + 1
            }
            #3) Update N.survive
            model$N.survive.F <<- model$N.F[2:n.year] - model$N.recruit.F #survivors are guys alive in year g-1 minus recruits in this year g
            model$calculate(ER.F.nodes) #update ER when N updated
          }
          #get proposed logprobs for N and y
          lp.proposed.N.M <- model$calculate(N.M.nodes[1])
          lp.proposed.N.F <- model$calculate(N.F.nodes[1])
          lp.proposed.N.recruit.M <- model$calculate(N.recruit.M.nodes)
          lp.proposed.N.recruit.F <- model$calculate(N.recruit.F.nodes)
          lp.proposed.y <- model$calculate(y.nodes[pick.idx]) #will always be 0

          lp.initial.total <- lp.initial.N.M + lp.initial.N.F + lp.initial.y + lp.initial.N.recruit.M + lp.initial.N.recruit.F
          lp.proposed.total <- lp.proposed.N.M + lp.proposed.N.F + lp.proposed.y + lp.proposed.N.recruit.M + lp.proposed.N.recruit.F

          #MH step
          log_MH_ratio <- lp.proposed.total - lp.initial.total
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["z.super",1] <<- model[["z.super"]]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            if(model$sex[pick]==0){
              mvSaved["N.M",1] <<- model[["N.M"]]
              mvSaved["N.survive.M",1] <<- model[["N.survive.M"]]
              mvSaved["N.recruit.M",1] <<- model[["N.recruit.M"]]
            }else{
              mvSaved["N.F",1] <<- model[["N.F"]]
              mvSaved["N.survive.F",1] <<- model[["N.survive.F"]]
              mvSaved["N.recruit.F",1] <<- model[["N.recruit.F"]]
            }
            mvSaved["ER.M",1] <<- model[["ER.M"]]
            mvSaved["ER.F",1] <<- model[["ER.F"]]
          }else{
            model[["z.super"]] <<- mvSaved["z.super",1]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            if(model$sex[pick]==0){
              model[["N.M"]] <<- mvSaved["N.M",1]
              model[["N.survive.M"]] <<- mvSaved["N.survive.M",1]
              model[["N.recruit.M"]] <<- mvSaved["N.recruit.M",1]
            }else{
              model[["N.F"]] <<- mvSaved["N.F",1]
              model[["N.survive.F"]] <<- mvSaved["N.survive.F",1]
              model[["N.recruit.F"]] <<- mvSaved["N.recruit.F",1]
            }
            model[["ER.M"]] <<- mvSaved["ER.M",1]
            model[["ER.F"]] <<- mvSaved["ER.F",1]
            #set these logProbs back
            model$calculate(N.M.nodes[1])
            model$calculate(N.F.nodes[1])
            model$calculate(N.recruit.M.nodes)
            model$calculate(N.recruit.F.nodes)
            model$calculate(y.nodes[pick.idx])
          }
        }
      }
    }
    
    #4) Finally, detected guy unobserved sex update
    for(i in 1:length(sex.up)){
      if(z.obs[sex.up[i]]==1){ #only do detected guys here
        i.idx <- seq(sex.up[i],M*n.year,M) #used to reference correct y nodes
        i.idx2 <- seq(sex.up[i],M*(n.year-1),M) #used to reference correct phi nodes 
        
        #get initial logProbs
        lp.initial.N.M <- model$getLogProb(N.M.nodes[1])
        lp.initial.N.F <- model$getLogProb(N.F.nodes[1])
        lp.initial.N.recruit.M <- model$getLogProb(N.recruit.M.nodes)
        lp.initial.N.recruit.F <- model$getLogProb(N.recruit.F.nodes)
        lp.initial.z <- model$getLogProb(z.nodes[sex.up[i]])
        lp.initial.y <- model$getLogProb(y.nodes[i.idx])
        #update N variables
        if(model$sex[sex.up[i]]==0){ #initial male
          #move this guy from N.M to N.F
          model$N.M <<- model$N.M - model$z[sex.up[i],]
          model$N.F <<- model$N.F + model$z[sex.up[i],]
          #move male recruit to female recruit
          if(model$z.start[sex.up[i]]>1){ #otherwise, this is year 1, nothing to change
            model$N.recruit.M[model$z.start[sex.up[i]]-1] <<- model$N.recruit.M[model$z.start[sex.up[i]]-1] - 1
            model$N.recruit.F[model$z.start[sex.up[i]]-1] <<- model$N.recruit.F[model$z.start[sex.up[i]]-1] + 1
          }
          # #update male, female survivors
          model$N.survive.M <<- model$N.M[2:n.year] - model$N.recruit.M
          model$N.survive.F <<- model$N.F[2:n.year] - model$N.recruit.F
        }else{ #initial female
          #move this guy from N.F to N.M
          model$N.F <<- model$N.F - model$z[sex.up[i],]
          model$N.M <<- model$N.M + model$z[sex.up[i],]
          #move female recruit to male recruit
          if(model$z.start[sex.up[i]]>1){ #otherwise, this is year 1, nothing to change
            model$N.recruit.F[model$z.start[sex.up[i]]-1] <<- model$N.recruit.F[model$z.start[sex.up[i]]-1] - 1
            model$N.recruit.M[model$z.start[sex.up[i]]-1] <<- model$N.recruit.M[model$z.start[sex.up[i]]-1] + 1
          }
          # #update male, female survivors
          model$N.survive.F <<- model$N.F[2:n.year] - model$N.recruit.F
          model$N.survive.M <<- model$N.M[2:n.year] - model$N.recruit.M
        }
        #update sex
        model$sex[sex.up[i]] <<- 1 - model$sex[sex.up[i]]
        #update phi nodes when sex changes
        model$calculate(phi.nodes[i.idx2])
        #get proposed logProbs
        lp.proposed.N.M <- model$calculate(N.M.nodes[1])
        lp.proposed.N.F <- model$calculate(N.F.nodes[1])
        lp.proposed.N.recruit.M <- model$calculate(N.recruit.M.nodes)
        lp.proposed.N.recruit.F <- model$calculate(N.recruit.F.nodes)
        lp.proposed.z <- model$calculate(z.nodes[sex.up[i]])
        lp.proposed.y <- model$calculate(y.nodes[i.idx])
        
        lp.initial.total <- lp.initial.z + lp.initial.y + lp.initial.N.recruit.M + lp.initial.N.recruit.F + lp.initial.N.M + lp.initial.N.F 
        lp.proposed.total <- lp.proposed.z + lp.proposed.y + lp.proposed.N.recruit.M + lp.proposed.N.recruit.F + lp.proposed.N.M + lp.proposed.N.F
        
        #MH step
        log_MH_ratio <- lp.proposed.total - lp.initial.total
        accept <- decide(log_MH_ratio)
        
        if(accept) {
          mvSaved["N.M",1] <<- model[["N.M"]]
          mvSaved["N.F",1] <<- model[["N.F"]]
          mvSaved["N.survive.M",1] <<- model[["N.survive.M"]]
          mvSaved["N.survive.F",1] <<- model[["N.survive.F"]]
          mvSaved["N.recruit.M",1] <<- model[["N.recruit.M"]]
          mvSaved["N.recruit.F",1] <<- model[["N.recruit.F"]]
          mvSaved["sex",1][sex.up[i]] <<- model[["sex"]][sex.up[i]]
          for(g in 1:(n.year-1)){
            mvSaved["phi",1][sex.up[i],g] <<- model[["phi"]][sex.up[i],g]
          }
        }else{
          model[["N.M"]] <<- mvSaved["N.M",1]
          model[["N.F"]] <<- mvSaved["N.F",1]
          model[["N.survive.M"]] <<- mvSaved["N.survive.M",1]
          model[["N.survive.F"]] <<- mvSaved["N.survive.F",1]
          model[["N.recruit.M"]] <<- mvSaved["N.recruit.M",1]
          model[["N.recruit.F"]] <<- mvSaved["N.recruit.F",1]
          model[["sex"]][sex.up[i]] <<- mvSaved["sex",1][sex.up[i]]
          for(g in 1:(n.year-1)){
            model[["phi"]][sex.up[i],g] <<- mvSaved["phi",1][sex.up[i],g]
          }
          model$calculate(y.nodes[i.idx])
          model$calculate(phi.nodes[i.idx2])
          model$calculate(z.nodes[sex.up[i]])
          model$calculate(N.recruit.M.nodes)
          model$calculate(N.recruit.F.nodes)
          model$calculate(N.M.nodes[1])
          model$calculate(N.F.nodes[1])
        }
      }
    }

    #copy back to mySaved to update logProbs.
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)

