rm(list = ls())
library(lavaan)
library(MASS)
library(MIIVsem)

N <- 1000
mu_latent <- c(0,0) #latent factor means #k=#of latent factors
varcov_latent <- matrix(c(1, .3, #var-cov matrix for latent factors
                          .3, 1), nrow = 2, byrow = T)
Lambda <- matrix(c(1, .8, .7, .7, 0, 0, 0, 0, #loadings
                   0, 0, 0, 0, 1, .8, .7, .7), nrow = 8, byrow = F) #p*k matrix.p=#of observed variables.
mu_epsilon <- rep(0, 8)
varcov_epsilon <- diag(1, nrow = 8)
var_epsilon_vector <- c(1-1^2, #standardized
                        1-.8^2,
                        1-.7^2,
                        1-.7^2,
                        1-1^2,
                        1-.8^2,
                        1-.7^2,
                        1-.7^2)
varcov_epsilon <- diag(c(1-1^2, #standardized
                         1-.8^2,
                         1-.7^2,
                         1-.6^2,
                         1-1^2,
                         1-.8^2,
                         1-.7^2,
                         1-.6^2),
                       nrow = 8)

set.seed(1234)
zeta <- mvrnorm(n = N,
                mu = mu_latent,
                Sigma = varcov_latent) #n*p matrix
#generate residuals
set.seed(1234)
epsilon <- mvrnorm(n = N,
                   mu =mu_epsilon,
                   Sigma = varcov_epsilon) #p*n matrix

#### Generate the observed variable series ####
y <- matrix(0, nrow = N, ncol = 8)
for (p in 1:nrow(y)){
  y[p, ] <- Lambda %*% zeta[p, ] + sqrt(var_epsilon_vector)*epsilon[p, ]
}
colnames(y) <- paste0('x', c(1:8))
head(y)


summary(miive(model = 'f1=~x1+x2+x3+x4
            f2=~x5+x6+x7+x8',
              data = y,
              var.cov = F))
