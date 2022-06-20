rm(list = ls())
library(lavaan)
library(MASS)
library(MIIVsem)

####simple generation###
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
              data = gendata), rsquare = T)

#####simulate 2 extra columns
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
                    0,.65,
                    .05,.05,
                    .05,.05,
                    .05,.05), nrow = 11, byrow = T)
varcov_epsilon <- diag(.25, nrow = 11)
seed <- 1
#data generation
set.seed(1234+seed)
zeta <- mvrnorm(n = N,
                mu = mean_latent,
                Sigma = varcov_latent) #generate latent factors

set.seed(1234+seed)
epsilon <- mvrnorm(n = N,
                   mu = rep(0, 11),
                   Sigma = varcov_epsilon) #n*p matrix, generate error terms

gendata <- matrix(0, nrow = 500, ncol = 11)
for (p in 1:nrow(gendata)){
  gendata[p, ] <- loading %*% zeta[p, ] + epsilon[p, ]
}
head(gendata)
#or
gendata <- t(loading %*% t(zeta)) + epsilon
head(gendata)
colnames(gendata) <- paste0('x', c(1:11))

heatmap(cor(gendata[,1:8]))

summary(miive(model = 'f1=~x1+x2+x3+x4
            f2=~x5+x6+x7+x8',
              data = gendata,
              var.cov = T))
cov(gendata[,1:8])

