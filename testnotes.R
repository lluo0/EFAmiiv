#######3.26.22########

#######Simulate 3########
model_simdata3 <- 'f1 =~ 1*z1 + .61 * z2 + .61*z3 + .61*z4
        f2=~ 1*z5 + .61*z6 + .61*z7 + .61*z8
        f1 ~~ .25*f2
z1 ~~ .3*z5
z2 ~~ .3*z6
z3 ~~ .3*z7
z4 ~~ .3*z8'
set.seed(1234)
data_simdata3 <- simulateData(model_simdata3, sample.nobs = 500)

EFAmiiv(data_simdata3, .05,'sargan+factorloading_R2',
          correlatedErrors = 'z1 ~~z5
          z2~~z6
          z3~~z7
          z4~~z8') #no error msg
EFAmiiv(data_simdata3, .05,'sargan+factorloading_R2') #ERROR msg

data <- data_simdata3
sigLevel <- .05
scalingCrit <- 'sargan+factorloading_R2'
correlatedErrors <- NULL

#error from miivsem shows: Error in processData(data, sample.cov, sample.mean, sample.nobs, ordered,  :
#object 'sample.sscp' not found
miive(model = 'f1=~z1+z3+z2+z5+z6\nf2=~z8+z3+z2+z5+z6\nf3=~z7+z3+z2+z5+z6\nf4=~z4+z3+z2+z5+z6',
      data = data_simdata3, var.cov = T)





#######empirical 1 #######
empirical1 <- lavaan::PoliticalDemocracy[,1:8]
head(empirical1)

EFAmiiv_cbnd(data = empirical1, sigLevel = .05, scalingCrit = 'sargan+factorloading_R2',
          correlatedErrors = 'y1 ~~y5
          y2~~y6
          y3~~y7
          y4~~y8')

EFAmiiv_cbnd(model = 'f1=~y1+y2+y3+y4
             f2=~y5+y6+y7+y8',
             data = empirical1, sigLevel = .05, scalingCrit = 'sargan+factorloading_R2',
             correlatedErrors = 'y1 ~~y5
          y2~~y6
          y3~~y7
          y4~~y8')

EFAmiiv_cbnd(model = 'f1=~y1+y2+y3+y4
             f2=~y5+y6+y7+y8',
             data = empirical1, sigLevel = .01, scalingCrit = 'sargan+factorloading_R2',
             correlatedErrors = 'y1 ~~y5
          y2~~y6
          y3~~y7
          y4~~y8')

EFAmiiv_cbnd(model = 'f1=~y1+y2+y3+y4+y5+y6+y7+y8',
             data = empirical1, sigLevel = .05, scalingCrit = 'sargan+factorloading_R2',
             correlatedErrors = 'y1 ~~y5
          y2~~y6
          y3~~y7
          y4~~y8')
#########empirical 2#######
#accdata <- read.table("/Users/lanluo/Downloads/access_raw.txt",header=F,sep=",")
accdata <- read.table("/Volumes/GoogleDrive/My Drive/EFAmiive paper/access_raw.txt",header=F,sep=",")
colnames(accdata) <- c(sapply(c(1:6), function(x) paste0('access', x)), sapply(c(1:6), function(x) paste0('easy', x)))

EFAmiiv_cbnd(data = accdata, sigLevel = .05, scalingCrit = 'sargan+factorloading_R2',
          correlatedErrors = 'access1~~easy1
      access2~~easy2
      access3~~easy3
      access4~~easy4
      access5~~easy5
      access6~~easy6')

EFAmiiv_cbnd(data = accdata, sigLevel = .01, scalingCrit = 'sargan+factorloading_R2',
             correlatedErrors = 'access1~~easy1
      access2~~easy2
      access3~~easy3
      access4~~easy4
      access5~~easy5
      access6~~easy6') #change the significance level to .01 leads to a clean 2 factor model

EFAmiiv_cbnd(data = accdata, sigLevel = .05, scalingCrit = 'sargan+factorloading_R2') #!underidentified model

data <- accdata
sigLevel <- .05
scalingCrit <- 'sargan+factorloading_R2'
##problem: atfer using var.cov = F, only prints out a partial model.

#########empirical 3#######
