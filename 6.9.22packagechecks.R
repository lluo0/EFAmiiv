
library(lavaan)
library(MIIVsem)
library(mvtnorm)
#####load functions####
#load functions in file:
#1, allfunctions_EFAmiiv.R
#2, allfunctions_EFAmiiv_multi.R
#3, function_getparam.R
#4, function_FAsim.R

#### simulation samples ####
#### sim1: simple crossload case ####
sim1 <- lapply(0:29, function(x) FAsim(n = 500,
                                       model = 'f1=~1*x1+.8*x2+.7*x3+.7*x4+.4*x7
            f2=~1*x5+.8*x6+.7*x7+.7*x8
                 f1 ~~ .6*f2',
                                      mean_latent = 0,
                                      seed = x, #30 reps
                                      var_obs = 1.25))

View(sim1)
EFAmiiv(data = sim1[[1]],
        sigLevel = .05,
        scalingCrit = 'sargan+factorloading_R2',
        correlatedErrors = NULL) #recovers true DGM
EFAmiiv(data = sim1[[2]],
        sigLevel = .05,
        scalingCrit = 'sargan+factorloading_R2',
        correlatedErrors = NULL) #does not crossload x7 because significant sargan even after crossloading, see below
miive('f1=~x1+x2+x3+x4+x7
      f2=~x5+x6+x7+x8', data = sim1[[2]], var.cov = T)

#### sim2: correlated errors between variables on top of crossloading ####
sim2 <- lapply(0:29, function(x) FAsim(n = 500,
                                       model = 'f1=~1*x1+.8*x2+.7*x3+.7*x4+.4*x7
            f2=~1*x5+.8*x6+.7*x7+.7*x8
                 f1 ~~ .6*f2
                                       x3 ~~ .4*x4',
                                       mean_latent = 0,
                                       seed = x, #30 reps
                                       var_obs = 1.25))

EFAmiiv(data = sim2[[1]],
        sigLevel = .05,
        scalingCrit = 'sargan+factorloading_R2',
        correlatedErrors = NULL) #creates additional factor for x3 and x4
EFAmiiv(data = sim2[[1]],
        sigLevel = .05,
        scalingCrit = 'sargan+factorloading_R2',
        correlatedErrors = 'x3 ~~ x4') #recovers the true DGM after including correlated errors

EFAmiiv_multi(model = 'f1=~x1+x2+x3+x4
              f2=~x5+x6+x7+x8', 
              data = sim2[[1]],
              sigLevel = .05,
              scalingCrit = 'sargan+factorloading_R2',
              correlatedErrors = NULL) #still creates additional factor for x3 and x4

EFAmiiv_multi(model = 'f1=~x1+x2+x3+x4
              f2=~x5+x6+x7+x8', 
              data = sim2[[1]],
              sigLevel = .05,
              scalingCrit = 'sargan+factorloading_R2',
              correlatedErrors = 'x3~~x4') #similar to when we start with no initial model

s1 <- EFAmiiv_multi_s1(model = 'f1=~x1+x2+x3+x4
              f2=~x5+x6+x7+x8', 
                       data = sim2[[1]],
                       sigLevel = .05,
                       scalingCrit = 'sargan+factorloading_R2',
                       correlatedErrors = NULL)#similar to when we start with no initial model
#### sim 3: 3 crossloadings for 3 factor model####
sim3 <- lapply(0:29, function(x) FAsim(n = 1000, #increased number of measurements because more complicated model
                                       model = 'f1=~1*x1+.8*x2+.7*x3+.7*x4+.4*x7
            f2=~1*x5+.8*x6+.7*x7+.7*x8+.4*x11
            f3=~1*x9+.8*x10+.7*x11+.6*x12+.3*x3
                 f1 ~~ .6*f2
                                       f1~~.5*f3
                                       f2~~.5*f3',
                                       mean_latent = 0,
                                       seed = x, #30 reps
                                       var_obs = 1.25))
EFAmiiv(data = sim3[[1]],
        sigLevel = .05,
        scalingCrit = 'sargan+factorloading_R2',
        correlatedErrors = NULL) #identical to true DGM
EFAmiiv(data = sim3[[2]],
        sigLevel = .05,
        scalingCrit = 'sargan+factorloading_R2',
        correlatedErrors = NULL) #identical to true DGM


