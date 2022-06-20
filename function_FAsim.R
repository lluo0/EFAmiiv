# FAsim <- function(N,
#                   mean_latent = 0,
#                   varcov_latent,
#                   loading,
#                   varcov_epsilon,
#                   seed = 0){
#   #setup
#   num_latent <- dim(loading)[2] # number of latent factors
#   num_measured <- dim(loading)[1] # number of observed variables
#
#   #data generation
#   set.seed(1234+seed)
#   if(num_latent==1){
#     zeta <- rnorm(n = N, mean = mean_latent, sd = sqrt(as.vector(varcov_latent)))
#   }
#   if(num_latent!=1){
#     zeta <- mvtnorm::rmvnorm(n = N,
#                              mean = mean_latent,
#                              sigma = varcov_latent) #generate latent factors
#   }
#   set.seed(1234+seed)
#   epsilon <- mvtnorm::rmvnorm(n = N,
#                      mean = rep(0, num_measured),
#                      sigma = varcov_epsilon) #n*p matrix, generate error terms
#
#   gendata <- t(loading%*%t(zeta))+epsilon #compute observed variables
#   if(!is.null(rownames(loading))){
#     colnames(gendata) <- rownames(loading)
#   }
#   return(gendata)
# }

####FAsim when specify a model directly#####
FAsim <- function(n,
                  model,
                  mean_latent = 0,
                  seed = 0,
                  var_obs = 1.25, #observed variables' error variance set to be 1.25 by default
                  var_error = NULL,
                  R2 = NULL){
  #get var-cov matrices for loadings, latent factors and errors
  modelparam <- getparam(model,
                         var_obs = var_obs, #observed variables' error variance set to be 1.25 by default
                         var_error = var_error,
                         R2 = R2)

  loading <- modelparam$Lambda
  varcov_latent <- modelparam$varcov_latent
  varcov_epsilon <- modelparam$varcov_epsilon

  num_latent <- dim(loading)[2] # number of latent factors
  num_measured <- dim(loading)[1] # number of observed variables

  #data generation
  set.seed(1234+seed)
  if(num_latent==1){
    zeta <- rnorm(n = n, mean = mean_latent, sd = sqrt(as.vector(varcov_latent)))
  }
  if(num_latent!=1){
    if(length(mean_latent)==1){
      zeta <- mvtnorm::rmvnorm(n = n,
                               mean = rep(mean_latent, num_latent),
                               sigma = varcov_latent) #generate latent factors
    }
    if(length(mean_latent)!=1){
      zeta <- mvtnorm::rmvnorm(n = n,
                               mean = mean_latent,
                               sigma = varcov_latent) #generate latent factors
    }
  }
  set.seed(1234+seed)
  epsilon <- mvtnorm::rmvnorm(n = n,
                              mean = rep(0, num_measured),
                              sigma = varcov_epsilon) #n*p matrix, generate error terms

  gendata <- t(loading%*%t(zeta))+epsilon #compute observed variables
  if(!is.null(rownames(loading))){
    colnames(gendata) <- rownames(loading)
  }
  return(as.data.frame(gendata))
}
#
# #####examples#####
#
# ########examples######
# #1
# model <- 'f1=~1*x1 + .8*x2 + .7*x3 + .7*x4 + .3*x7
# f2=~ 1*x5 + .8*x6 + .7*x7 + .7*x8
# f1 ~~ .5*f2
# f1 ~~ .4*f3
# f2 ~~ .6*f3'
# testparam <- getparam(model, R2 = c(.8, .7, .7, .7,
#                                       .6, .8, .7, .7,
#                                       .8, .7, .7, .7))
# # testsim <- FAsim(N = 200,
# #                  mean_latent = rep(0, 3),
# #                  varcov_latent =  testparam$varcov_latent,
# #                  loading = testparam$Lambda,
# #                  varcov_epsilon = testparam$varcov_epsilon,
# #                  seed = 1)
#
# summary(cfa('f1=~x1+x2+x3+x4
#             f2=~x5+x6+x7+x8
#             f3=~x9+x10+x11+x12', testsim), rsquare =T)
# #2
# model <- 'f1=~1*x1 + .8*x2 + .7*x3 + .7*x4 + .3*x7
# f2=~ 1*x5 + .8*x6 + .7*x7 + .7*x8
# f1 ~~ .5*f2'
# testparam <- getparam(model, R2 = rep(.8, 8))
# testsim <- FAsim(N = 1000,
#                  mean_latent = rep(0, 2),
#                  varcov_latent =  testparam$varcov_latent,
#                  loading = testparam$Lambda,
#                  varcov_epsilon = testparam$varcov_epsilon,
#                  seed = 1)
#
# summary(cfa('f1=~x1+x2+x3+x4
#             f2=~x5+x6+x7+x8', testsim), rsquare =T)
# summary(miive('f1=~x1+x2+x3+x4+x7
#             f2=~x5+x6+x7+x8', testsim, var.cov = T))
