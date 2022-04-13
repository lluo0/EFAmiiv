##tests and debugs##

########3.29.2022##########

##
#m3: 2 factor model w/ x7 crossload on f1
sm3 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .5*x7
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
        f1 ~~ .5*f2'
sim3 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim3[[p]] <- simulateData(sm3, sample.nobs = 1000)
}

EFAmiiv(data = sim3[[1]],
        sigLevel = .05,
        scalingCrit = 'sargan+factorloading_R2',
        correlatedErrors = NULL)

EFAmiiv_multi(model = 'f1=~x1+x2+x3+x4+x5+x6+x7+x8',
              data = sim3[[1]],
              sigLevel = .05,
              scalingCrit = 'sargan+factorloading_R2',
              correlatedErrors = NULL)

EFAmiiv_multi(model='f1+~x1+x2+x3+x4+x5+x6+x7+x8',
              data=sim3[[1]],
              sigLevel = 0.05,
              scalingCrit = 'sargan+factorloading_R2',
              correlatedErrors=NULL)

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
EFAmiiv(data = sim10[[1]],
        sigLevel = .05,
        scalingCrit = 'sargan+factorloading_R2',
        correlatedErrors = 'x2~~x3')
