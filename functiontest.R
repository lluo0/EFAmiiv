#tests
#load packages
library(MIIVsem)
library(psych)
library(lavaan)
##run all the functions in these files: allfunctions_EFAmiiv.R, allfunction_EFAmiiv_multi.R, and EFAmiiv_combined.R

##then test on simulation data in the simdata.R file - need to run the functions to generate the simulated data first.
##for each simulation there's multiple replications, use simdataname[[n]] to use the nth individual simulated dataset.
  ##for example, for simulation 3, run the following code
#m3: 2 factor model w/ x7 crossload on f1
sm3 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .5*x7
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
        f1 ~~ .5*f2'
sim3 <- list()
for (p in 1:30){ #so there's 30 reps
  set.seed(123.4+p)
  sim3[[p]] <- simulateData(sm3, sample.nobs = 1000)
}
#then use sim3[[n]] to test the function's performance:
EFAmiiv(data = sim3[[1]], #no need to manually test all 30reps, just choose a few
        sigLevel = .05, #feel free to change the significance level threshold here, usually we use .05, sometimes can change to .01 for empirical data
        scalingCrit = 'sargan+factorloading_R2',
  #there are multiple other options for the scalingCrit argument
  #see the scaling indicator selection section in the allfunctions.EFAmiiv.R file.
  #Feel free to play around with this argument and change it to other value and see if there's any error message.
        correlatedErrors = NULL) #can leave this as NULL for most cases

##please also test the performance/any error messages when including an initial model.
  #for example, for simulation 3:
EFAmiiv(model = 'f1=~x1+x2+x3+x4
        f2=~x5+x6+x7+x8', #this is almost identical to the true data generation model, except we omit the x7 on f1 - this is to test if the function can recover that
        #note that for model specification, no parameter coefficient included.
        data = sim3[[1]],
        sigLevel = .05,
        scalingCrit = 'sargan+factorloading_R2',
        correlatedErrors = NULL)
  #feel free to play around the initial model specification - it doesn't have to be identical to the true data generation model
  #but the dimensionality, aka the number of factors, should be less than the true model.
  #for example, you could start with a one factor model like:
EFAmiiv(model = 'f1=~x1+x2+x3+x4+x5+x6+x7+x8',
        data = sim3[[1]],
        sigLevel = .05,
        scalingCrit = 'sargan+factorloading_R2',
        correlatedErrors = NULL)

##for simulations with correlated errors in the data generation model, please test the performance when both including/excluding them.
  #for example, for simulation 10:
#m10: 2 factor, with a crossloading and a pair of correlated error
sm10 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .6*x7
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
      f1 ~~ .5*f2
        x2 ~~ .3*x3'
sim10<- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim10[[p]] <- simulateData(sm10, sample.nobs = 1000)
}

EFAmiiv(data = sim3[[1]],
        sigLevel = .05,
        scalingCrit = 'sargan+factorloading_R2',
        correlatedErrors = 'x2~~x3') #include x2~~x3 or leave it as NULL. note no coefficient here like in the model specification.
