rm(list = ls())
library(lavaan)
library(MASS)
library(MIIVsem)

describe(sim3[[1]])

sm3 <- 'f1 =~ 1*x1 + .8 * x2 + .7*x3 + .7*x4 + .5*x7
        f2=~ 1*x5 + .7*x6 + .6*x7 + .6*x8
        f1 ~~ .5*f2
        x1 ~~ (1-1^2)*x1
        x2 ~~ (1-.8^2)*x2
        x3 ~~ (1-.7^2)*x3
        x4 ~~ (1-.7^2)*x4
        x5 ~~ (1-1^2)*x5
        x6 ~~ (1-.7^2)*x6
        x7 ~~ (1-.5^2-.6^2)*x7
        x8 ~~ (1-.6^2)*x8
'
sim3 <- list()
for (p in 1:30){
  set.seed(123.4+p)
  sim3[[p]] <- simulateData(sm3, sample.nobs = 1000)
}
describe(sim3[[1]])

summary(sem(model = 'f1=~x1+x2+x3+x4+x7
            f2=~x5+x6+x7+x8', data = sim3[[1]]))

N <- 1000
mu_latent <- c(0,0) #latent factor means #k=#of latent factors
varcov_latent <- matrix(c(1, .3, #var-cov matrix for latent factors
                .3, 1), nrow = 2, byrow = T)
Lambda <- matrix(c(1, .8, .7, .7, 0, 0, 0, 0, #loadings
                   0, 0, 0, 0, 1, .8, .7, .7), nrow = 8, byrow = F) #p*k matrix.p=#of observed variables.
mu_epsilon <- rep(0, 8)
varcov_epsilon <- diag(.01, nrow = 8)
varcov_epsilon <- diag(c(1-1^2, #standardized
                         1-.8^2,
                         1-.7^2,
                         1-.7^2,
                         1-1^2,
                         1-.8^2,
                         1-.7^2,
                         1-.7^2),
                       nrow = 8)
#generate latent factors
library(MASS)
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
  y[p, ] <- Lambda %*% zeta[p, ] + epsilon[p, ]
}
colnames(y) <- paste0('x', c(1:8))
head(y)
#identical to the following
gendata <- t(Lambda%*%t(zeta) + t(epsilon)) #transpose of p*n matrix
dim(gendata)
cov(gendata)
gendata <- cbind(gendata, 1:N)
colnames(gendata) <- c(paste0('x', c(1:8)), 'Time')
head(gendata)
gendata_reshape <- reshape2::melt(as.data.frame(gendata), id = "Time")

ggplot(gendata_reshape, aes(Time, value)) + geom_line() +
  theme_bw() + ylab("") +
  facet_wrap(~ variable, ncol = 3) +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


EFAmiiv(data = as.data.frame(gendata))

summary(cfa(model = 'f1=~x1+x2+x3+x4
            f2=~x5+x6+x7+x8',
            data = gendata))
summary(miive(model = 'f1=~x1+x2+x3+x4+x7
            f2=~x5+x6+x7+x8',
              data = gendata,
              var.cov = T))


simEFA <- function(N, alpha, Phi, Lambda, Theta){

}


