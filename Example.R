#----------------------------------------------------------------------------------------------------#######
#----------- This is an example of estimating the most cost-effective individulized -----------------#######  
#----------- treatment rule using the package "OptimalCEITR"                        -----------------#######
#----------------------------------------------------------------------------------------------------#######

library(survival);library(rpart);library(dplyr);library(cubature)
setwd("C:/Users/cryst/Box Sync/Crystal's Thesis Shared Folder/Manuscript1-OptimalCEITR/FinalRevisionFolder")
source("./OptimalCEITR.R")

##----------------------------- Simulate covariates, treatment, and outcome NMB --------------------------##
N      <- 100
tau    <- 3         # restriction time point: 3 year
lambda <- 100000    # willingness-to-pay threshold

set.seed(99)
# Baseline covariates
x1 <- rnorm(N,0,1)
x2 <- rnorm(N,1,1)
x3 <- rnorm(N,0,1)

# Treatment 
u0 <- 0.1*x2
p0 <- 1/(1+exp(-u0))
A  <- rbinom(N,1,p0)

# Survival time (ST)
Main  <- as.matrix(data.frame(x1,x2,x3,A))
betaF <- c(0.1,0.3,0.3,0.2)
gammaF<- 0.6
ST    <- exp(Main%*%betaF+gammaF*x2*A)
C     <- runif(N,0,15)
FT    <- pmin(ST,C)
event <- ifelse(ST<=C,1,0)
RFT   <- pmin(FT,rep(tau,N))

# Cumulative cost (CC)
gammaC <- 0.8
beta   <- 1000*exp(Main%*%betaF+gammaC*x2*A)
YC     <- rgamma(N,shape=2.5,scale=beta)
CC     <- YC*RFT

# CE outcome
Y <- lambda*RFT-CC;summary(Y)


##---------------------------------------------- Analysis --------------------------------------------------##
# estimated outcomes
dat <- data.frame(x1,x2,x3,A,event,FT,RFT,CC,Y)
muT <- Reg.mu.surv(ST=FT, event=event, A=A, Xs=data.frame(x1,x3), Ms=x2, tau=tau, data=dat)$mus.survRT
muC <- Reg.mu.norm(CC=CC, A=A, Xs=data.frame(x1,x3), Ms=x2, data=dat)$mus.norm

# estimated propensity scores
est_ps <- EstPS(A=A,X=cbind(x1,x2,x3))

# estimated censoring weights
complete <- ifelse((dat$event==1 | dat$FT>=tau),1,0)   # complete:having event|contribute full follow-up of interest
cen.wm   <- glm(complete ~ A*x1+A*x2+A*x3, family = "binomial")
cen.w    <- predict(cen.wm, type = "response")

# estimated AIPW classification weights
ITR.aipw <- ITR_AIPW(Q = RFT, C = CC, A = A, est_ps = est_ps, cen = event, est_cen = cen.w,
                     lambda = lambda, muQ.reg = muT, muC.reg = muC)
aipw.Y <- ITR.aipw$delta.Y
aipw.ITRCE <- ITR.aipw$ITR_aipwCE

# estimated individualized treatment rule 
dat <- data.frame(x1,x2,x3,aipw.Y,aipw.ITRCE)
comp.dat <- dat[complete.cases(dat)==T,]
fit.aipw <- rpart(aipw.ITRCE ~ x1+x2+x3, weights=comp.dat$aipw.Y, method="class",data=comp.dat)
comp.dat$g.aipw <- as.numeric(predict(fit.aipw, type="class"))-1;table(comp.dat$g.aipw)
imp.g.aipw <- fit.aipw$variable.importance;imp.g.aipw


