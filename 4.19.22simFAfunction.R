rm(list = ls())
library(lavaan)
library(MASS)
library(MIIVsem)

#solution: generate 2 more extra variables to discard!
########one factor model######
L11 <- 1
L21 <- .8
L31 <- .75
L41 <- .7
L51 <- .65
L61 <- .6
L71 <- .55
L81 <- .3
L91 <- .05
L101 <- .05
# L12 <- 1
# L22 <- .8
# L32 <- .7
# L42 <- .7
# L52 <- .6
# L62 <- .6
# L72 <- .5
# L82 <- .5

N <- 1000
set.seed(1234)
zeta <- rnorm(n=N, mean = 0, sd = 1)
Lambda <- c(L11, L21, L31, L41, L51, L61, L71, L81, L91,L101)

v1 <- 1-L11^2
v2 <- 1-L21^2
v3 <- 1-L31^2
v4 <- 1-L41^2
v5 <- 1-L51^2
v6 <- 1-L61^2
v7 <- 1-L71^2
v8 <- 1-L81^2
v9 <- 1-L91^2
v10 <- 1-L101^2
mu_epsilon <- rep(0,10)
varcov_epsilon <- diag(c(v1, v2, v3, v4, v5, v6, v7, v8,v9,v10), nrow = 10)
set.seed(1234)
epsilon <- mvrnorm(n = N,
                   mu =mu_epsilon,
                   Sigma = varcov_epsilon) #n*p matrix

sim_1fac <- t(Lambda%*%t(zeta))+epsilon
corrplot::corrplot(cor(sim_1fac), method = 'number')

head(sim_1fac)

pairs(sim_1fac[,c(1,8)])
cor(sim_1fac)


########two factor model#######

L11 <- 1
L21 <- .8
L31 <- .7
L41 <- .6
L51 <- 0
L61 <- 0
L71 <- .3
L81 <- 0
L91 <- .05
L101 <- .05

L12 <- 0
L22 <- 0
L32 <- 0
L42 <- 0
L52 <- 1
L62 <- .8
L72 <- .7
L82 <- .6
L92 <- .05
L102 <- .05

N <- 1000
set.seed(12341)
varcov_zeta <- matrix(c(1,.3,
                        .3,1), nrow = 2, byrow = F)
zeta <- mvrnorm(n = N,
                mu = c(0,0),
                Sigma = varcov_zeta)
Lambda <- cbind(c(L11, L21, L31, L41, L51, L61, L71, L81, L91,L101),
                c(L12, L22, L32, L42, L52, L62, L72, L82, L92,L102))

v1 <- 1-L11^2-L12^2
v2 <- 1-L21^2-L22^2
v3 <- 1-L31^2-L32^2
v4 <- 1-L41^2-L42^2
v5 <- 1-L51^2-L52^2
v6 <- 1-L61^2-L62^2
v7 <- 1-L71^2-L72^2
v8 <- 1-L81^2-L82^2
v9 <- 1-L91^2-L92^2
v10 <- 1-L101^2-L102^2
mu_epsilon <- rep(0,10)
varcov_epsilon <- diag(c(v1, v2, v3, v4, v5, v6, v7, v8,v9,v10), nrow = 10)
set.seed(12341)
epsilon <- mvrnorm(n = N,
                   mu =mu_epsilon,
                   Sigma = varcov_epsilon) #n*p matrix

sim_2facCL <- t(Lambda%*%t(zeta))+epsilon
corrplot::corrplot(cor(sim_2facCL), method = 'circle')

colnames(sim_2facCL) <- c(paste0('x', c(1:10)))

summary(cfa(model = 'f1=~x1+x2+x3+x4
            f2=~x5+x6+x7+x8',
            data = sim_2facCL))
summary(miive(model = 'f1=~x1+x2+x3+x4+x7
            f2=~x5+x6+x7+x8',
              data = sim_2facCL,
              var.cov = T))
##simple 2 factor model
testsim <- FAsim(N = 500,
      mean_latent = c(0,0),
      varcov_latent = matrix(c(1,.4,
                               .4,1), nrow = 2, byrow = T),
      loading = matrix(c(1,0,
                         .8,0,
                         .7,0,
                         .6,0,
                         0,1,
                         0,.8,
                         0,.7,
                         .4,.6), nrow = 8,byrow = T),
      varcov_epsilon = diag(.25, nrow = 8),
      seed = 1
      )
colnames(testsim) <- paste0('x', c(1:8))
EFAmiiv(as.data.frame(testsim))

testparam <- getparam(model = 'f1=~1*x1 + .8*x2 + .7*x3 + .6*x4
         f2=~1*x5 + .8*x6 + .7*x7 + .6*x8
         f1~~ .5*f2', epsilon_var = .25, latent_var = 1)

N <- 500
mean_latent = c(0,0)
varcov_latent <- testparam$varcov_latent
loading <- testparam$Lambda
varcov_epsilon <- testparam$varcov_epsilon
seed = 1
testsim <- FAsim(N = 500,
                 mean_latent = c(0,0),
                 varcov_latent = testparam$varcov_latent,
                 loading = testparam$Lambda,
                 varcov_epsilon = testparam$varcov_epsilon,
                 seed = 1)

#data generation
set.seed(1234+seed)
zeta <- mvrnorm(n = N,
                mu = mean_latent,
                Sigma = varcov_latent) #generate latent factors

set.seed(1234+seed)
epsilon <- mvrnorm(n = N,
                   mu = rep(0, 8),
                   Sigma = varcov_epsilon) #n*p matrix, generate error terms

gendata <- matrix(0, nrow = 500, ncol = 8)
for (p in 1:nrow(gendata)){
  gendata[p, ] <- loading %*% zeta[p, ] + epsilon[p, ]
}
colnames(gendata) <- rownames(loading)

heatmap(cor(testsim))
summary(miive(model = 'f1=~x1+x2+x3+x4
            f2=~x5+x6+x7+x8',
              data = testsim,
              var.cov = T))

EFAmiiv(as.data.frame(gendata))


FAsim <- function(N,
                  mean_latent,
                  varcov_latent,
                  loading,
                  varcov_epsilon,
                  seed = 0){
  #setup
  num_latent <- dim(loading)[2] # number of latent factors
  num_measured <- dim(loading)[1] # number of observed variables
  # num_measured_extra <- max(num_measured/2, 5)
  # extraloading <- matrix(.05, nrow = num_measured_extra, ncol = num_latent)
  # loading_add <- rbind(loading, extraloading) #add 5 more additional variables to be discarded later
  # #setup var-cov matrix for epsilons
  # var_epsilon <- 1-rowSums(loading_add^2)
  # varcov_epsilon <- varcov_epsilon
  # varcov_epsilon_add <- cbind(rbind(varcov_epsilon,
  #                                   matrix(0, nrow = num_measured_extra, ncol = nrow(varcov_epsilon))),
  #                             matrix(0, nrow = num_measured_extra+num_measured, ncol = num_measured_extra))
  # varcov_epsilon <- matrix(0,
  #                          nrow = num_measured+num_measured_extra,
  #                          ncol = num_measured+num_measured_extra)
  # diag(varcov_epsilon) <- var_epsilon #need to add a stage to update based on correlated errors
  # diag(varcov_epsilon_add)[(num_measured+1):(num_measured+num_measured_extra)] <- var_epsilon[(num_measured+1):(num_measured+num_measured_extra)]

  #data generation
  set.seed(1234+seed)
  zeta <- mvrnorm(n = N,
                  mu = mean_latent,
                  Sigma = varcov_latent) #generate latent factors

  set.seed(1234+seed)
  epsilon <- mvrnorm(n = N,
                     mu = rep(0, num_measured),
                    # mu = rep(0,num_measured+num_measured_extra),
                     #Sigma = varcov_epsilon_add
                    Sigma = varcov_epsilon) #n*p matrix, generate error terms

  # gendata <- t(loading_add%*%t(zeta))+epsilon #compute observed variables
  # gendata <- gendata[,c(1:num_measured)] #get rid of extra generations
  gendata <- t(loading%*%t(zeta))+epsilon #compute observed variables
  colnames(gendata) <- rownames(loading)
  return(gendata)
}
corrplot::corrplot(cor(gendata), method = 'number')

corrplot::corrplot(cor(gendata[,c(1:num_measured)]), method = 'number')
heatmap(cov(gendata))
colnames(gendata) <- paste0('x', c(1:8))
summary(cfa(model = 'f1=~x1+x2+x3+x4
            f2=~x5+x6+x7+x8',
            data = gendata))
summary(miive(model = 'f1=~x1+x2+x3+x4+x7
            f2=~x5+x6+x7+x8',
              data = gendata,
              var.cov = T))

summary(miive(model = 'f1=~x1+x2+x3+x4+x7
            f2=~x5+x6+x7+x8',
              data = gendata,
              var.cov = T))

EFAmiiv(as.data.frame(gendata))

summary(miive(model = 'f1=~x1+x2+x3+x4
            f2=~x5+x6+x7+x8',
              data = gendata,
              var.cov = T))

#functions to extract specific parameters from model specification

#models would be latent =~ , latent ~~ latent, and observed ~~ observed

getparam <- function(model, epsilon_var, latent_var=1){
  latentFac <- observedVar <- vector()
  ######setup#####
  modelstring <- strsplit(model, "\n")[[1]]
  #get the number of equations specifying loadings. no correlated errors here.
  latentEq <- as.vector(which(sapply(modelstring, function(u) grepl('=~', u))))
  #get the names of the latent factors. !use gsub to get rid of spaces
  latentFac <- gsub(' ','',as.vector(sapply(modelstring[latentEq], function(i) strsplit(i, '=~')[[1]][1])))
  # #the number of equations with correlated errors for latent factors
  # corrEq_latentFac <- as.vector(which(sapply(modelstring,
  #                                            function(u) grepl(paste(latentFac, collapse = '|'), u))))[-latentEq]
  # #the number of equations with correlated errors for observed variables
  # corrEq_observedVar <- c(1:length(modelstring))[-c(latentEq, corrEq_latentFac)]

  #####get names and loadings for observed variables#####
  observedVar_list <- lapply(modelstring[latentEq],
                              function(u) gsub(".*[=~]([^.]+)", "\\1", u)) #get strings after =~
  observedVar_list <- lapply(observedVar_list,
                        function(u) strsplit(u, '+', fixed = T)[[1]]) #separate strings by +
  # observedVar <- lapply(observedVar_list,
  #                       function(u) as.vector(sapply(u,
  #                                          function(i) gsub(".*[*]([^.]+)", "\\1", i))))
  observedVar <- lapply(observedVar_list,
                        function(u) as.vector(sapply(u,
                         function(i) gsub(' ', '', sapply(i, #get rid of spaces
                           function(m) gsub(".*[*]([^.]+)", "\\1", m)))))) #get strings after *
  observedVar_loadings <- lapply(observedVar_list,
                                 function(u) as.numeric(sapply(u,
                                  function(i) gsub(' ', '', sapply(i, #get rid of spaces
                                   function(m) gsub("[*].*$", "", m)))))) #get strings before *
  observedVar_unique <- unique(unlist(observedVar))
  ######create the loading matrix#####
  Lambda <- matrix(0,
                   nrow = length(observedVar_unique),
                   ncol = length(latentFac))
  rownames(Lambda) <- observedVar_unique
  colnames(Lambda) <- latentFac
  for(n in 1:ncol(Lambda)){
    for(p in 1:nrow(Lambda)){
      if(rownames(Lambda)[p] %in% observedVar[[n]]){
        ordernum <- which(observedVar[[n]] == rownames(Lambda)[p]) #get the corresponding loadings' order in that list
        Lambda[p,n] <- as.numeric(observedVar_loadings[[n]][ordernum])
      }
    }
  }

  #####get the equation numbers for var-cov matrices#####
  #the number of equations with correlated errors for observed variables
  corrEq_observedVar <- as.vector(which(sapply(modelstring,
                                               function(u) grepl(paste(observedVar_unique, collapse = '|'), u))))[-latentEq]
  #the number of equations with correlated errors for latent factors
  corrEq_latentFac <- c(1:length(modelstring))[-c(latentEq, corrEq_observedVar)]

  ######get covariance of latent factors######
  if(length(corrEq_latentFac)==0){ #no correlation between latent factors
    varcov_latent <- matrix(0, nrow = length(latentFac), ncol = length(latentFac))
    colnames(varcov_latent) <- rownames(varcov_latent) <- latentFac
  }
  if(length(corrEq_latentFac!=0)){ #get the correlation between latent factors based on specified model
    latent_list <- lapply(modelstring[corrEq_latentFac],
                          function(u) strsplit(u, '~~', fixed = T)[[1]]) #separate strings by +
    for(n in 1:length(latent_list)){
      names(latent_list)[n] <- gsub(".*[*]([^.]+)", "\\1", latent_list[[n]][2]) #name each list the other latent factor
      latent_list[[n]][1] <- gsub(' ', '', latent_list[[n]][1]) #get rid of spaces
      latent_list[[n]][2] <- as.numeric(gsub("[*].*$", "", latent_list[[n]][2])) #leave the 2nd element the covariance
    }
    ######create latent factor var-cov matrix#####
    varcov_latent <- diag(0, nrow = length(latentFac))
    colnames(varcov_latent) <- rownames(varcov_latent) <- latentFac
    for(n in 1:length(latent_list)){
      for(p in 1:nrow(varcov_latent)){
        if(names(latent_list)[n]==rownames(varcov_latent)[p]){
          for(i in 1:ncol(varcov_latent)){
            # varcov_latent[p,i] <- varcov_latent[i,p] <- ifelse(colnames(varcov_latent)[i]==latent_list[[n]][1],
            #                              as.numeric(latent_list[[n]][2]), 0)
            varcov_latent[i,p] <-  ifelse(colnames(varcov_latent)[i]==latent_list[[n]][1],
                                          as.numeric(latent_list[[n]][2]), varcov_latent[i,p])
            varcov_latent[p,i] <-  ifelse(colnames(varcov_latent)[i]==latent_list[[n]][1],
                                          as.numeric(latent_list[[n]][2]), varcov_latent[p,i])
          }
        }
      }
    }

  }
  diag(varcov_latent) <- latent_var #variance of latent factors set to be 1 if not specified otherwise

  ######get covariance of observed variables errors#####
  if(length(corrEq_observedVar)==0){
    varcov_epsilon <- matrix(0, nrow = length(observedVar_unique), ncol = length(observedVar_unique) )
    colnames(varcov_epsilon) <- rownames(varcov_epsilon) <- observedVar_unique
  }
  if(length(corrEq_observedVar)!=0){
    epsilon_list <- lapply(modelstring[corrEq_observedVar],
                           function(u) strsplit(u, '~~', fixed = T)[[1]]) #separate strings by +
    for(n in 1:length(epsilon_list)){
      names(epsilon_list)[n] <- gsub(' ','',gsub(".*[*]([^.]+)", "\\1", epsilon_list[[n]][2])) #name each list the other latent factor
      epsilon_list[[n]][1] <- gsub(' ', '', epsilon_list[[n]][1]) #get rid of spaces
      epsilon_list[[n]][2] <- as.numeric(gsub("[*].*$", "", epsilon_list[[n]][2])) #leave the 2nd element the covariance
    }
    ######create latent factor var-cov matrix#####
    varcov_epsilon <- diag(0, nrow = length(observedVar_unique))
    colnames(varcov_epsilon) <- rownames(varcov_epsilon) <- observedVar_unique
    for(n in 1:length(epsilon_list)){
      for(p in 1:nrow(varcov_epsilon)){
        if(names(epsilon_list)[n]==rownames(varcov_epsilon)[p]){
          for(i in 1:ncol(varcov_epsilon)){
            varcov_epsilon[p,i] <- ifelse(colnames(varcov_epsilon)[i]==epsilon_list[[n]][1],
                                          as.numeric(epsilon_list[[n]][2]), varcov_epsilon[p,i])
            varcov_epsilon[i,p] <- ifelse(colnames(varcov_epsilon)[i]==epsilon_list[[n]][1],
                                          as.numeric(epsilon_list[[n]][2]), varcov_epsilon[i,p])
          }
        }
      }
    }
  }

  # diag(varcov_epsilon) <- 1-rowSums(Lambda^2) #standardize
  diag(varcov_epsilon) <- epsilon_var
  finalobj <- list(Lambda = Lambda,
                   varcov_latent = varcov_latent,
                   varcov_epsilon = varcov_epsilon)
  return(finalobj)
}

getparam('factor1 =~ 1*m190+.8*m290+.7*m390+.6*x490+.2*m795
         factor2=~1*m595+.8*m695+.7*m795+.8*m895
         factor3=~1*m990+.8*m1090+.7*m1190+.6*m1290
         factor1~~.3*factor2
         factor2~~.4*factor3
         factor1 ~~ .3*factor3
         m1190~~.2*m1290')

getparam(model = 'f1=~1*x1 + .8*x2 + .7*x3 + .6*x4 + .4*x8
         f2=~1*x5 + .8*x6 + .7*x7 + .6*x8
         f1~~ .5*f2', epsilon_var = .25, latent_var = 1)

model <- 'factor1 =~ 1*m190+.8*m290+.7*m390+.6*m490+.2*m795
         factor2=~1*m595+.8*m695+.7*m795+.8*m895
         factor3=~1*m990+.8*m1090+.7*m1190+.6*m1290
         factor1~~.3*factor2
         factor2~~.4*factor3
         factor1 ~~ .3*factor3
         m1190~~.2*m1290'


model <- 'f1=~1*x1 + .8*x2 + .7*x3 + .7*x4 + .3*x7
f2=~ 1*x5 + .8*x6 + .7*x7 + .7*x8
f3=~ 1*x9 + .8*x10 + .7*x11 + .6*x12
f1 ~~ .3*f2
f1 ~~ .4*f3
f2 ~~ .35*f3
x2 ~~ .3*x3'

epislon_var <- .25

testparam <- getparam(model)

testdata <- FAsim(N = 1000,
                  mean_latent = rep(0,3),
                  varcov_latent = testparam$varcov_latent,
                  loading = testparam$Lambda,
                  varcov_epsilon = testparam$varcov_epsilon,
                  seed = 0)
corrplot::corrplot(cor(testdata), method = 'number')

cov(testdata)

EFAmiiv(as.data.frame(testdata), correlatedErrors = 'x2~~x3')
s1 <- step1_EFAmiiv(as.data.frame(testdata), .05, correlatedErrors = 'x2~~x3')
s2 <- step2_EFAmiiv(s1, as.data.frame(testdata), .05, 'sargan+factorloading_R2')
s3 <- stepN_EFAmiiv(s2, as.data.frame(testdata), .05, 'sargan+factorloading_R2')
s4 <- stepN_EFAmiiv(s3, as.data.frame(testdata), .05, 'sargan+factorloading_R2')

stepPrev <- s2
