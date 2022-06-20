rm(list = ls())
library(lavaan)
library(MASS)
library(MIIVsem)

#### original 2 factor model that had issues with covariances ####
N <- 500
mean_latent <- c(0,0)
varcov_latent <- matrix(c(1,.5,
                          .5, 1), nrow = 2)
loading <- matrix(c(1,0,
                    .8,0,
                    .7,0,
                    .6,0,
                    0,1,
                    0,.85,
                    0,.75,
                    0,.65), nrow = 8, byrow = T)
varcov_epsilon <- diag(.25, nrow = 8)
seed <- 1
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
head(gendata)
#or
gendata <- t(loading %*% t(zeta)) + epsilon
head(gendata)
colnames(gendata) <- paste0('x', c(1:8))

heatmap(cor(gendata))
cov(gendata)

summary(miive(model = 'f1=~x1+x2+x3+x4
            f2=~x5+x6+x7+x8',
              data = gendata,
              var.cov = T))
summary(cfa(model = 'f1=~x1+x2+x3+x4
            f2=~x5+x6+x7+x8',
            data = gendata), rsquare = T) #r2 seems fine, neither too large or too small.
# note that r2 for x8 should be roughly identical to x4, but x8's r2 is obviously too large based on the true DGM.

#### simple 1 factor model ####
N <- 500
mean_latent <- 0
varcov_latent <- 1
loading <- matrix(c(1,
                    .8,
                    .7,
                    .6), nrow = 4, byrow = T)
varcov_epsilon <- diag(c(.25, .16, .1225, .09), nrow = 4)
seed <- 0
#data generation
set.seed(1234+seed)
zeta <- mvrnorm(n = N,
                mu = mean_latent,
                Sigma = varcov_latent) #generate latent factors

set.seed(1234+seed)
epsilon <- mvrnorm(n = N,
                   mu = rep(0, 4),
                   Sigma = varcov_epsilon) #n*p matrix, generate error terms
#compute observed variables
gendata <- t(loading %*% t(zeta)) + epsilon
head(gendata)
colnames(gendata) <- paste0('x', c(1:4))
cov(gendata)
#the covariance matrix of the observed data should be: Lambda*phi*t(Lambda) + theta
loading%*%varcov_latent%*%t(loading) + varcov_epsilon

heatmap(cor(gendata))


summary(cfa('f1=~x1+x2+x3+x4', data = gendata), rsquare = T)
summary(miive('f1=~x1+x2+x3+x4', data = gendata, var.cov = T))



#generate data based on the covariance matrix of observed data
y <- mvrnorm(n = N,
             mu = rep(0, 4),
             Sigma = (loading%*%varcov_latent%*%t(loading) + varcov_epsilon))
cov(y)
colnames(y) <- paste0('x', c(1:4))
summary(cfa('f1=~x1+x2+x3+x4', data = y), rsquare = T)
summary(miive('f1=~x1+x2+x3+x4', data = y, var.cov = T))


#data generation
set.seed(1234+seed)
# zeta <- mvtnorm::rmvnorm(n = N,
#                 mean = c(mean_latent),
#                 sigma = varcov_latent) #generate latent factors
zeta <- rnorm(n = N, mean = 0, sd = sqrt(varcov_latent))

set.seed(1234+seed)
epsilon <- mvtnorm::rmvnorm(n = N,
                   mean = rep(0, 4),
                   sigma = varcov_epsilon) #n*p matrix, generate error terms
#compute observed variables
gendata <- t(loading %*% t(zeta)) + epsilon
head(gendata)
colnames(gendata) <- paste0('x', c(1:4))
cov(gendata)

#####2#####
N <- 500
mean_latent <- c(0,0)
varcov_latent <- matrix(c(1,.5,
                          .5, 1), nrow = 2)
loading <- matrix(c(1,0,
                    .8,0,
                    .7,0,
                    .6,0,
                    0,1,
                    0,.85,
                    0,.75,
                    0,.65), nrow = 8, byrow = T)
varcov_epsilon <- diag(.25, nrow = 8)

set.seed(1234+seed)
zeta <- mvtnorm::rmvnorm(n = N, mean = mean_latent, sigma = varcov_latent)
cov(zeta)
set.seed(1234+seed)
epsilon <- mvtnorm::rmvnorm(n = N,
                            mean = rep(0, 8),
                            sigma = varcov_epsilon) #n*p matrix, generate error terms
cov(epsilon)
#compute observed variables
gendata <- t(loading %*% t(zeta)) + epsilon

cov(gendata)
loading%*%varcov_latent%*%t(loading) + varcov_epsilon

colnames(gendata) <- paste0('x', c(1:8))

summary(miive(model = 'f1=~x1+x2+x3+x4
            f2=~x5+x6+x7+x8',
              data = gendata,
              var.cov = T))
summary(cfa(model = 'f1=~x1+x2+x3+x4
            f2=~x5+x6+x7+x8',
            data = gendata), rsquare = T)
