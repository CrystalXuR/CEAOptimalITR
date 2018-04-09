#### This simulation data is built based on the description in papaer titled ###########
#### "Estimating optimal treatment regimen from a classification perspective" ##########



################################## Data Simulation ################################
### Scenario I & II: g.opt functions are the only difference 
NDataSets = 1000
# Different sample sizes for each data set
N1 = 200
N2 = 500
N3 = 1000
set.seed(1)
DataI.200 = list(NA,200)
DataII.200 = list(NA,200)
for (i in 1:NDataSets){
  N = N1
  # 5 covariates ~ random normal
  x1 = rnorm(N,0.3,0.5)
  x2 = rnorm(N,0.3,0.8)
  x3 = rnorm(N,1,2)
  x4 = rnorm(N,2,5)
  x5 = rnorm(N,0.4,1)
  
  # binary treatment 
  reg = -0.1+0.5*x1+0.5*x2
  u = exp(reg)/(1+exp(reg))
  A = rbinom(N,1,u)
  
  # outcome Y with epsilang ~ N(0, 1)
  g.optI = I(x1>-0.545)*I(x2<0.545)
  YI = exp(2+0.25*x1+0.25*x2-0.25*x5-0.5*(A-g.optI)^2)+rnorm(N,0,1)
  DataI.200[[i]] = data.frame(cbind(x1,x2,x3,x4,x5,A,g.optI,YI))
  
  g.optII = I(x1>x2)
  YII = exp(2+0.25*x1+0.25*x2-0.25*x5-0.5*(A-g.optII)^2)+rnorm(N,0,1)
  DataII.200[[i]] = data.frame(cbind(x1,x2,x3,x4,x5,A,g.optII,YII))
}


DataI.500 = list(NA,500)
DataII.500 = list(NA,500)
for (i in 1:NDataSets){
  N = N2
  # 5 covariates ~ random normal
  x1 = rnorm(N,0.3,0.5)
  x2 = rnorm(N,0.3,0.8)
  x3 = rnorm(N,1,2)
  x4 = rnorm(N,2,5)
  x5 = rnorm(N,0.4,1)
  
  # binary treatment 
  reg = -0.1+0.5*x1+0.5*x2
  u = exp(reg)/(1+exp(reg))
  A = rbinom(N,1,u)
  
  # outcome Y with epsilang ~ N(0, 1)
  g.optI = I(x1>-0.545)*I(x2<0.545)
  YI = exp(2+0.25*x1+0.25*x2-0.25*x5-0.5*(A-g.optI)^2)+rnorm(N,0,1)
  DataI.500[[i]] = data.frame(cbind(x1,x2,x3,x4,x5,A,g.optI,YI))
  
  g.optII = I(x1>x2)
  YII = exp(2+0.25*x1+0.25*x2-0.25*x5-0.5*(A-g.optII)^2)+rnorm(N,0,1)
  DataII.500[[i]] = data.frame(cbind(x1,x2,x3,x4,x5,A,g.optII,YII))
}


DataI.1000 = list(NA,1000)
DataII.1000 = list(NA,1000)
for (i in 1:NDataSets){
  N = N3
  # 5 covariates ~ random normal
  x1 = rnorm(N,0.3,0.5)
  x2 = rnorm(N,0.3,0.8)
  x3 = rnorm(N,1,2)
  x4 = rnorm(N,2,5)
  x5 = rnorm(N,0.4,1)
  
  # binary treatment 
  reg = -0.1+0.5*x1+0.5*x2
  u = exp(reg)/(1+exp(reg))
  A = rbinom(N,1,u)
  
  # outcome Y with epsilang ~ N(0, 1)
  g.optI = I(x1>-0.545)*I(x2<0.545)
  YI = exp(2+0.25*x1+0.25*x2-0.25*x5-0.5*(A-g.optI)^2)+rnorm(N,0,1)
  DataI.1000[[i]] = data.frame(cbind(x1,x2,x3,x4,x5,A,g.optI,YI))
  
  g.optII = I(x1>x2)
  YII = exp(2+0.25*x1+0.25*x2-0.25*x5-0.5*(A-g.optII)^2)+rnorm(N,0,1)
  DataII.1000[[i]] = data.frame(cbind(x1,x2,x3,x4,x5,A,g.optII,YII))
}


################################# Data Implementation #############################
#### Contrast function, C.hat: three choices: regression, IPWE, AIPWE
## Truth using counterfactual data
test.data = data.frame(DataI.1000[[1]]); names(test.data);dim(test.data) #200 8
u1 = exp(2+0.25*test.data$x1+0.25*test.data$x2-0.25*test.data$x5-0.5*(1-test.data$g.optI)^2)
u0 = exp(2+0.25*test.data$x1+0.25*test.data$x2-0.25*test.data$x5-0.5*(0-test.data$g.optI)^2)
# True E{Y*(g.opt)} = 8.26
mean(u1*test.data$g.optI+u0*(1-test.data$g.optI)) 


############ CART to choose the optimal treatment regime ########################
### Estimate contrast function using regression - incorrect model specification
reg.fit = lm(YI~A*x1+A*x2+A*x3+A*x4+A*x5, data = test.data); summary(reg.fit)  

N = N3
A = rep(1,N)
treated.data = cbind(test.data[, c(1:5,7,8)], A)
A = rep(0,N)
Nontreated.data = cbind(test.data[, c(1:5,7,8)], A)

reg.out1 = predict(reg.fit, type = "response",newdata = treated.data)
reg.out0 = predict(reg.fit, type = "response", newdata = Nontreated.data)
C.reg = reg.out1-reg.out0
Z.hat = ifelse(C.reg>0,1,0)
W.hat = abs(C.reg)
class.data = data.frame(cbind(test.data,Z.hat,W.hat)); dim(class.data) #200 10

# Weighted classification using CART
fit.reg<-rpart(Z.hat ~ x1+x2+x3+x4+x5, weights=W.hat, method="class",data=class.data)

# optimal treatment from regression 
gopt.reg = as.numeric(predict(fit.reg, type = 'class', data = class.data))
gopt.reg = ifelse(gopt.reg==2,1,0)

# Estimate the optimal outcome
A = gopt.reg
reg.data = cbind(test.data[, c(1:5,7,8)], A)
#reg.opt1 = predict(reg.fit, type = "response",newdata = reg.data)
reg.opt1 = exp(2+0.25*test.data$x1+0.25*test.data$x2-0.25*test.data$x5-0.5*(1-gopt.reg)^2)
mean(reg.opt1) #6.63







### Using AIPW to estimate the contrast function
# Correct ps
test.data$A = as.factor(test.data$A)
PS.model = glm(A~x1+x2, family = binomial, data = test.data)
PS.1 = predict(PS.model, type = "response")
PS.0 = 1-PS.1

# Incorrect ps
PS.model2 = glm(A~1, family = binomial, data = test.data)
InPS.1 = predict(PS.model2, type = "response")
InPS.0 = 1-InPS.1

### estimate conditional means by regression
REG2<-Reg.mu(Y=test.data$YI,As=test.data$A,H=cbind(test.data$x1,test.data$x2,test.data$x3,test.data$x4,test.data$x5))
mus2.reg<-REG2$mus.reg
fit = REG2$RegModel

### calculate AIPW adaptive contrasts C and working orders l (regimen)
# a1 is for ACWL-C1 and a2 is for ACWL-C2
# input outcome and treatment vectors, estimated propensity and conditional means
CLs2.a<-CL.AIPW(Y=test.data$YI,A=test.data$A,pis.hat=cbind(PS.1, PS.0),mus.reg=mus2.reg)
C2.a1<-CLs2.a$C.a1
#l2.a<-CLs2.a$l.a            #optimal treatment regime

# when PS is incorreclty specified
CLs2.a<-CL.AIPW(Y=test.data$YI,A=test.data$A,pis.hat=cbind(InPS.1, InPS.0),mus.reg=mus2.reg)
InC2.a1<-CLs2.a$C.a1


### Weighted classification using CART
fit.AIPW<-rpart(Z.hat ~ x1+x2+x3+x4+x5, weights=C2.a1, method="class",data=test.data)
gopt.AIPW = as.numeric(predict(fit.AIPW,type = "class", data=test.data))
gopt.AIPW = ifelse(gopt.AIPW==2,1,0)

# when PS is incorreclty specified
Infit.AIPW<-rpart(Z.hat ~ x1+x2+x3+x4+x5, weights=InC2.a1, method="class",data=test.data)
Ingopt.AIPW = as.numeric(predict(Infit.AIPW,type = "class", data=test.data))
Ingopt.AIPW = ifelse(Ingopt.AIPW==2,1,0)


### Estimate the optimal outcome
A = gopt.AIPW
AIPW.data = cbind(test.data[, c(1:5,7,8)], A)
#AIPW.opt1 = predict(reg.fit, type = "response",newdata = AIPW.data)
AIPW.opt1 = exp(2+0.25*test.data$x1+0.25*test.data$x2-0.25*test.data$x5-0.5*(1-gopt.AIPW)^2)
mean(AIPW.opt1) #6.67

# when PS is incorreclty specified
A = Ingopt.AIPW
InAIPW.data = cbind(test.data[, c(1:5,7,8)], A)
#InAIPW.opt1 = predict(reg.fit, type = "response",newdata = InAIPW.data)
InAIPW.opt1 = exp(2+0.25*test.data$x1+0.25*test.data$x2-0.25*test.data$x5-0.5*(1-Ingopt.AIPW)^2)
mean(InAIPW.opt1) #6.64

# test the optimal regimen aganist the truth 
test1 = ifelse(test.data$g.optI==gopt.reg, 1,0);table(test1)
test2 = ifelse(test.data$g.optI==gopt.AIPW, 1,0);table(test2)# 912(1000)
test3 = ifelse(test.data$g.optI==Ingopt.AIPW, 1,0);table(test3) # 903(1000) 
## In both AIPW cases, the outcome model is incorrectly specified, 
## but when the PS is also incorrect specified, the 
## estimate optimal regimen has a higher misclassification rate.










#----------------------------------------------------------------#
#------------------ optimalclass code ---------------------------#
#----------------------------------------------------------------#

## install required packages  ##
install.packages("R.utils");library(R.utils)
#install.packages("DynTxRegime")
library(DynTxRegime);library(rpart)

## Load the source code   ##
source("C:/Users/u0918455/Box Sync/A_Crystal's Thesis/1 - Optimization-Classification/2 - Simulations & implementation/Simulation Codes and Results/Optimal Class Code/optimalClasscx.r")


# consider 2 PS models: one is correct, one is false
moPropen <- buildModelObj(model = ~ test.data$x1+test.data$x2,
                              solver.method = 'glm',
                              solver.args = list('family'='binomial'),
                              predict.method = 'predict.glm',
                              predict.args = list(type='response'))

# Define the classification model.
moClass <- buildModelObj(model = ~ test.data$x1+test.data$x2+test.data$x5,
                             solver.method = 'rpart',
                             solver.args = list(method="class"),
                             predict.args = list(type='class'))

# Define the cenosring model.
moCensor <- buildModelObj(model = ~ test.data$x1+test.data$x2+test.data$A,
                          solver.method = 'glm',
                          solver.args = list('family'='binomial'),
                          predict.method = 'predict.glm',
                          predict.args = list(type='response'))

### IPWE estimator
# both PS and classification model are correct
test.data$CI = exp(2+2.5*x1+2.5*x2-0.25*x5-5*(A-g.optI)^2)+rnorm(N,0,1)
censored = rbinom(1000, 1, 0.5)

estIPWE.TT <- .optimalClasscx(moPropen = moPropen, moClass = moClass,moCensor = moCensor, 
                                  lambda = 50000,censored = censored, 
                                  data = test.data, responseE = test.data$YI,responseC = test.data$CI, txName = "A",iter = 0L)







# # PS is correct but outcome model is wrong
# estIPWE.FT <- optimalClass(moPropen = moPropenFalse, moClass = moClassTrue,
#                         data = test.data, response = test.data$YI, txName = "A",iter = 0L)
# 
# # PS is wrong but outcome model is correct
# estIPWE.TF <- optimalClass(moPropen = moPropenTrue, moClass = moClassFalse,
#                         data = test.data, response = test.data$YI, txName = "A",iter = 0L)
# 
# # both PS and outcome are wrong
# estIPWE.FF <- optimalClass(moPropen = moPropenFalse, moClass = moClassFalse,
#                         data = test.data, response = test.data$YI, txName = "A",iter = 0L)

# Output(optimal regimen) from IPW
IPW.Regime.TT = as.numeric(optTx(estIPWE.TT, newdata = test.data))
# IPW.Regime.FT = as.numeric(optTx(estIPWE.FT, newdata = test.data))
# IPW.Regime.TF = as.numeric(optTx(estIPWE.TF, newdata = test.data))
# IPW.Regime.FF = as.numeric(optTx(estIPWE.FF, newdata = test.data))

# ## Test the misclassification rate against the truth
# test4a = ifelse(test.data$g.optI==IPW.Regime.TT, 1,0);table(test4a) # 958(1000) ## Compare to the other code ('biom---'), OptimalClass perfoms better
# test4b = ifelse(test.data$g.optI==IPW.Regime.FT, 1,0);table(test4b) # 767
# test4c = ifelse(test.data$g.optI==IPW.Regime.TF, 1,0);table(test4c) # 928
# test4d = ifelse(test.data$g.optI==IPW.Regime.FF, 1,0);table(test4d) # 597
# 
# ## Estimate the optimal outcome 
# summary(estIPWE.TT)$value  #8.386
# summary(estIPWE.FT)$value  #7.683
# summary(estIPWE.TF)$value  #8.581
# summary(estIPWE.FF)$value  #8.247  # when both models are incorrect, the outcome even higher~ ~




# ### AIPWE estimator
# # Create modelObj object for main effect component
# #moMain <- buildModelObj(model = ~ test.data$A*test.data$x1+test.data$A*test.data$x2+test.data$A*test.data$x3+test.data$A*test.data$x4+test.data$A*test.data$x5, solver.method = 'lm')
# moMain <- buildModelObj(model = ~ test.data$A+test.data$x1+test.data$x2+test.data$x5, solver.method = 'lm')
# 
# # Create modelObj object for contrast component
# moCont <- buildModelObj(model = ~ test.data$x1+test.data$x2+test.data$x3+test.data$x4+test.data$x5, solver.method = 'lm')
# 
# estAIPWE <- optimalClass(moPropen = moPropenTrue, moMain = moMain,
#                          moCont = moCont, moClass = moClassTrue,
#                          data = test.data, response = test.data$YI, txName = "A",iter = 0L)
# 
# # output from AIPW
# summary(estAIPWE)$value
















