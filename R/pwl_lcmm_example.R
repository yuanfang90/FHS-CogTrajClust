rm(list=ls())
library(dplyr)
library(haven)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)

# require(devtools)
# install_version("lcmm", version = "1.9.3", repos = "http://cran.us.r-project.org")
library(lcmm)

source("stepwise.R")
source("gridsearch3.R")

set.seed(10, kind = "Mersenne-Twister", normal.kind = "Inversion",sample.kind = "Rejection")

model.data$age.1 <- model.data$age.centered
model.data$age.31 <- with(model.data,ifelse(age.1 > -0.5,age.1 + 0.5,0)) # 65
model.data$age.32 <- with(model.data,ifelse(age.1 > 0,age.1,0)) # 70
model.data$age.33 <- with(model.data,ifelse(age.1 > 0.5,age.1 - 0.5,0)) # 75
model.data$age.34 <- with(model.data,ifelse(age.1 > 1,age.1 - 1,0)) #80
model.data$age.35 <- with(model.data,ifelse(age.1 > 1.5,age.1 - 1.5,0)) #85
model.data$age.36 <- with(model.data,ifelse(age.1 > 2,age.1 - 2,0)) #90

####################################
####### 3 change points

#### G=1
set.seed(10, kind = "Mersenne-Twister", normal.kind = "Inversion",sample.kind = "Rejection")
ef.fact.m1.pwl.3 <- hlme(EF_EI ~ age.1+age.32+age.34+age.36+SEX+EDUcomp, random = ~ age.1+age.32+age.34+age.36, subject = "framid", verbose = TRUE, data = model.data)
minit <<- ef.fact.m1.pwl.3

#### G=2
set.seed(10, kind = "Mersenne-Twister", normal.kind = "Inversion",sample.kind = "Rejection")
m2.binit <- gridsearch3(nrep=30,niter=10,minit=ef.fact.m1.pwl.3,G=2,mod.dat=model.data)
ef.fact.m2.pwl.3 <- hlme(EF_EI ~ age.1+age.32+age.34+age.36+SEX+EDUcomp, random = ~ age.1+age.32+age.34+age.36, subject = "framid", data = model.data, ng = 2, verbose = TRUE, mixture = ~ age.1+age.32+age.34+age.36+SEX+EDUcomp,B=m2.binit)

#### G=3
set.seed(10, kind = "Mersenne-Twister", normal.kind = "Inversion",sample.kind = "Rejection")
m3.binit <- gridsearch3(nrep=30,niter=10,minit=ef.fact.m1.pwl.3,G=3,mod.dat=model.data)
ef.fact.m3.pwl.3 <- hlme(EF_EI ~ age.1+age.32+age.34+age.36+SEX+EDUcomp, random = ~ age.1+age.32+age.34+age.36, subject = "framid", data = model.data, ng = 3, verbose = TRUE, mixture = ~ age.1+age.32+age.34+age.36+SEX+EDUcomp,B=m3.binit)

#### G=4
set.seed(10, kind = "Mersenne-Twister", normal.kind = "Inversion",sample.kind = "Rejection")
m4.binit <- gridsearch3(nrep=30,niter=10,minit=ef.fact.m1.pwl.3,G=4,mod.dat=model.data)
ef.fact.m4.pwl.3 <- hlme(EF_EI ~ age.1+age.32+age.34+age.36+SEX+EDUcomp, random = ~ age.1+age.32+age.34+age.36, subject = "framid", data = model.data, ng = 4, verbose = TRUE, mixture = ~ age.1+age.32+age.34+age.36+SEX+EDUcomp,B=m4.binit)

########################################
### stepwise function
#### G=1
set.seed(10, kind = "Mersenne-Twister", normal.kind = "Inversion",sample.kind = "Rejection")
ef.fact.m1.step <- multcp.step(mod=ef.fact.m1.pwl.3,nrep=NA,niter=NA,n.cp.in=3)

#### G=2
set.seed(10, kind = "Mersenne-Twister", normal.kind = "Inversion",sample.kind = "Rejection")
ef.fact.m2.step <- multcp.step(mod=ef.fact.m2.pwl.3,nrep=30,niter=10,n.cp.in=3)

#### G=3
set.seed(10, kind = "Mersenne-Twister", normal.kind = "Inversion",sample.kind = "Rejection")
ef.fact.m3.step <- multcp.step(mod=ef.fact.m3.pwl.3,nrep=30,niter=10,n.cp.in=3)

#### G=4
set.seed(10, kind = "Mersenne-Twister", normal.kind = "Inversion",sample.kind = "Rejection")
ef.fact.m4.step <- multcp.step(mod=ef.fact.m4.pwl.3,nrep=30,niter=10,n.cp.in=3)

