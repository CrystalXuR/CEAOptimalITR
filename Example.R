#------------------------------------------------------------------------------------#######
#----------- This is an example of CEA using this package        --------------------#######
#----------- Outcome: treatment A                                --------------------#######
#----------- Predictor: covariate Xs                             --------------------#######
#------------------------------------------------------------------------------------#######
library(CEAOptimalITR); library(survival);library(rpart)
#------------------- Simulate covariates, treatment, and 2 outcomes -----------------#
N=100
tau=3         # restriction time point = 3 year
set.seed(99)
x1 = rnorm(N,0,1)
x2 = rnorm(N,1,1)
x3 = rnorm(N,0,1)

u0 = 0.1*x1
p0 = 1/(1+exp(-u0))
A = rbinom(N,1,p0)

# survival time (ST)
Main = as.matrix(data.frame(x1,x2,x3,A))
betaF=c(0.1,0.3,0.3,0.4); gammaF=0.3
ST = exp(Main%*%betaF+gammaF*x2*A);summary(ST)
RST = pmin(ST,rep(tau,N))   
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
lambda=100000
Y = lambda*FT-CC;summary(Y)


##-------------------- Analysis -------------------------##
#library(R.utils);library(survival);library(rpart)
#sourceDirectory("C:/Users/u0918455/Box/A_Crystal's Thesis/1 - Optimization-Classification/2 - Simulations & implementation/Simulation Codes and Results/1-CEanalysisITRcode/CE_ITR_SourceCode")

dat = data.frame(x1,x2,x3,A,FT,event,cen,CC,Y)
muT.reg.exp = Reg.mu.exp(ST=FT, A=A, Xs=as.matrix(data.frame(x1,x3)), Ms=x2, data=dat)
muTpre = muT.reg.exp
muT = matrix(NA,N,2)
for (j in 1:2){
  muT[,j] = pmin(rep(tau,N),muTpre[,j])
}

muC = Reg.mu.norm(CC=CC, A=A, Xs=as.matrix(data.frame(x1,x3)), Ms=x2, data=dat)

# estimated CE outcome
muY.reg.norm = lambda*muT-muC

# estimated PS
est_ps = EstPS(A=A,Xs=cbind(x1,x3))

# create censoring weight
cen.wm = glm(event ~ A+x1+x2+x3, family = "binomial");cen.w = predict(cen.wm, type = "response")

###------- calculate AIPW estimates -------###
ITR.aipw=ITR_AIPW(ST = FT, CC = CC, A = A, est_ps=est_ps, cen = event, est_cen = cen.w,lambda = lambda,
                  muT.reg = muT, muC.reg = muC)
aipw.Y=ITR.aipw$delta.Y
aipw.ITRCE=ITR.aipw$ITR_aipwCE

# AIPW contrasts based on CE
fit.aipw=rpart(aipw.ITRCE ~ x1+x2+x3, weights=aipw.Y, method="class",
               control = rpart.control(minsplit = 10, cp = 0.01,xval = 10))
g.aipw=as.numeric(predict(fit.aipw, type="class"))-1;table(g.aipw)
imp.g.aipw = fit.aipw$variable.importance;imp.g.aipw



