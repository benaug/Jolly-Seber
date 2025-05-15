# Jolly-Seber
Jolly-Seber MCMC samplers in Nimble


This repository contains a Jolly-Seber model and MCMC sampler that allows yearly per capita recruitment to be modeled as Poisson while allowing individual effects on survival, detection, or other things. This is achieved by representing survival and recruitment simultaneously in two sets of objects, the individual by year z matrix indicating years in the population, and N and N.recruit, the vectors for total abundance and numbers of recruits per year.

This approach has some efficiency benefits and perhaps costs (not sure, yet). In Nimble, it DEFINITELY has RAM benefits for data sets with massive numbers of z's. Using the typical approach, Nimble will assign a sampler to every single z and this can eat up a lot of RAM. In my approach, there is 1 sampler assigned to all z's. 

1) the dSurvival distribution only calculates the log-likelihood for survival for the relevant years (not before birth or the year after death year)

2) for individuals detected at least once, I update z.start and z.stop one at a time with a categorical sampler. The benefit of this approach is that mixing should be improved because we can turn on/off multiple z's for one individual per iteration. In fact, there could be cases where the typical approach updating one z at a time per individual cannot fully explore the posterior, at least, in theory. This new approach would allow you to "jump over" any severe bottlenecks between regions of the posterior. However, the categorical sampler may get very slow with many primary periods (e.g., years) due to the categorical sampler used. But it still may be more efficient than an "update one z at a time" approach. This should be investigated. An alternative proposal would be to propose an entire z vector at a time. This would avoid the slowness of the categorical update calculations, but I believe the acceptance rates for many individuals could be very low, at least using the same approach for proposing z vectors below. Here, I am thinking of individuals that recruited in later years--most of their proposals may be too early to be accepted.

3) for individuals never detected, whether in the super population or not: I propose entire z vectors for each individual. I propose the recruitment occasion, z.start, (including the possibility of being there in the first year). Then, I simulate survival from the proposed recruitment occasion. Then, I use a Metropolis-Hastings update. 

4) handling the superpopulation: I use count prior data augmentation https://github.com/benaug/SCR-Count-Prior-Data-Augmentation


Final Note: Ruth King suggested to me that the birth and death primary occasions could be updated instead of the individual z's right before I worked on this. I would like to think I would have thought of this myself, but can't be sure! So, I will credit Ruth for this observation. Of course, she bears no responsibility for the correctness of my full algorithm or code!
