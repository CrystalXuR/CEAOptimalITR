#------------------------------------------------------------------------------------#######
#----------- This is an example of CEA using this package        --------------------#######
#----------- Outcome: treatment A                                --------------------#######
#----------- Predictor: covariate Xs                             --------------------#######
#------------------------------------------------------------------------------------#######

#------------------- Simulate covariates, treatment, and 2 outcomes -----------------#
N=100
set.seed(99)
x1 = rnorm(N,0,1)
x2 = rnorm(N,0.5,1)
x3 = rnorm(N,0,1)

u0 = 0.1*x1
p0 = 1/(1+exp(-u0))
A = rbinom(N,1,p0)

# survival time (ST)
Main = as.matrix(data.frame(x1,x2,x3,A))
betaF=c(0.1,0.3,0.3,0.4); gammaF=0.3
ST = exp(Main%*%betaF+gammaF*x2*A);summary(ST)
RST = pmin(ST,rep(3,N))       # restriction time point = 3 year
C=runif(N,0,10)
FT = pmin(RST,C);summary(FT)
cen = ifelse(RST<=C,1,0)      # censoring indicator 
event = as.numeric(I(cen==1));table(event)

# Cumulative cost (CC)
gammaC=0.1
beta = 100*exp(Main%*%betaF+gammaC*x2*A)
YC = rgamma(N,shape=200,scale=beta);summary(YC)
CC = YC*FT;summary(CC)

# CE outcome
lambda=50000
Y = lambda*FT-CC;summary(Y)


##-------------------- Analysis -------------------------##
#library(R.utils);library(survival);library(rpart)
#sourceDirectory("C:/Users/cryst/Box Sync/A_Crystal's Thesis/1 - Optimization-Classification/2 - Simulations & implementation/Simulation Codes and Results/1-CEanalysisITRcode/CE_ITR_SourceCode")

dat = data.frame(x1,x2,x3,A,FT,event,cen,CC,Y)
muT.reg.exp = Reg.mu.exp(ST=FT, A=A, Xs=as.matrix(data.frame(x1,x3)), Ms=x2, data=dat)
muC.reg.norm = Reg.mu.norm(CC=CC, A=A, Xs=as.matrix(data.frame(x1,x3)), Ms=x2, data=dat)

# estimated CE outcome
muY.reg.norm = lambda*muT.reg.exp-muC.reg.norm

# estimated PS
est_ps = EstPS(A=A,Xs=cbind(x1,x3))

# create censoring weight 
cen.wm = glm(event ~ A, family = "binomial");cen.w = predict(cen.wm, type = "response") 

###------- calculate IPW estimates -------###
ITR.ipw=ITR_IPW(ST = FT, CC = CC, A = A, est_ps=est_ps, cen = event, est_cen = cen.w,lambda = lambda)
ipw.Y=ITR.ipw$delta.Y                     ## individual treatment effect on CE 
ipw.ITRCE=ITR.ipw$ITR_ipwCE               ## estimated regime based on CE

###------- calculate AIPW estimates -------###
ITR.aipw=ITR_AIPW(ST = FT, CC = CC, A = A, est_ps=est_ps, cen = event, est_cen = cen.w,lambda = lambda,
                  muT.reg = muT.reg.exp, muC.reg = muC.reg.norm)
aipw.Y=ITR.aipw$delta.Y
aipw.ITRCE=ITR.aipw$ITR_aipwCE

# IPW contrasts based on CE
fit.ipw=rpart(ipw.ITRCE ~ x1+x2+x3, weights=ipw.Y, method="class",
              control = rpart.control(minsplit = 10, cp = 0.01,xval = 10))
g.ipw=as.numeric(predict(fit.ipw, type="class"))-1;table(g.ipw)
imp.g.ipw = fit.ipw$variable.importance;imp.g.ipw

# AIPW contrasts based on CE
fit.aipw=rpart(aipw.ITRCE ~ x1+x2+x3, weights=aipw.Y, method="class",
               control = rpart.control(minsplit = 10, cp = 0.01,xval = 10))
g.aipw=as.numeric(predict(fit.aipw, type="class"))-1;table(g.aipw)
imp.g.aipw = fit.aipw$variable.importance;imp.g.aipw



