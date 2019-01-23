#####---------------------------------------------------------------------------#####
#####---------------- Functions for cost-effectiveness Analysis ----------------#####
#####---------------------------------------------------------------------------#####

# Estimate PS using logistic regression
# Outcome: treatment A 
# Predictor: covariate Xs
EstPS = function(A,Xs){
  if(ncol(as.matrix(A))!=1) stop("Multiple treatments are not allowed")
  if(length(A)!= nrow(as.matrix(Xs))) stop("Non-conformable between outcome and predictors")
  if(length(unique(A))<=1) stop("Outcome levels are less than 2")
  level.A = sort(unique(A))

  dat = data.frame(A,Xs)
  PSfit = glm(A ~., family="binomial", data=dat)
  p = predict(PSfit,dat,type = "response")
  PS = data.frame((1-p),p);names(PS) = paste("A=",level.A,sep="");PS  
}

# Estimate survival time using Cox PH model 
# Outcome: survival time ST (continuous)
# Event indicator: event
# Predictor: treatment A, main covariates Xs, and interaction Ms(*A).
Reg.mu.Cox = function(ST,A,Xs,Ms,event,data){
  #if(length(A) != dim(Xs)[1]|length(A) != dim(Ms)[1]|dim(Xs)[1] != dim(Ms)[1]) stop("Non-conformable between treatment and covariates")
  N = length(A)
  Xs = as.matrix(Xs)
  Ms = as.matrix(Ms)
  LA = length(unique(A))    
   if(LA<2) stop("Treatment levels are less than 2")

    cox.ph.model = coxph(formula = Surv(time = ST, event = event) ~ Xs+A*Ms, method="breslow",data=data)
    mus.CoxT = matrix(NA,N,LA)
    for(j in 1:LA){
      newdat = data.frame(Xs,A=rep(sort(unique(A))[j],N),Ms)
      survObj = survfit(cox.ph.model,newdata=newdat)
      surv.out = data.frame(summary(survObj)$table)
      mus.CoxT[,j] = surv.out$X.rmean     
    }
  
  outCoxT = list(mus.CoxT, cox.ph.model);names(outCoxT) = c("mus.CoxT","CoxPH");outCoxT
}

# Estimate cumulative cost using Lognormal model  
# Outcome: cumulative cost CC (continuous)
# Predictor: treatment A, main covariates Xs, and interaction Ms(*A).
Reg.mu.norm = function(CC,A,Xs,Ms,data){
  #if(length(A) != dim(Xs)[1]|length(A) != dim(Ms)[1]|dim(Xs)[1] != dim(Ms)[1]) stop("Non-conformable between treatment and covariates")
  N = nrow(as.matrix(A))
  Xs = as.matrix(Xs)
  Ms = as.matrix(Ms)
  LA = length(unique(A))    
   if(LA<2) stop("Treatment levels are less than 2")
  
    RegNormal = glm(CC ~ Xs+A*Ms, family = "gaussian"(link="log"),data = data)
    mus.norm = matrix(NA,N,LA)
    for(j in 1:LA){
      newdat = data.frame(Xs,A=rep(sort(unique(A))[j],N),Ms)
      mus.norm[,j] = predict(RegNormal,newdata=newdat)
    }

  out.norm = list(mus.norm, RegNormal);names(out.norm) = c("mus.norm","LogNormal");out.norm
}

# Estimate cumulative cost using Gamma (identity) model  
# Outcome: cumulative cost CC (continuous)
# Predictor: treatment A, main covariates Xs, and interaction Ms(*A).
Reg.mu.gam = function(CC,A,Xs,Ms,data){
  #if(length(A) != dim(Xs)[1]|length(A) != dim(Ms)[1]|dim(Xs)[1] != dim(Ms)[1]) stop("Non-conformable between treatment and covariates")
  N = nrow(as.matrix(A))
  Xs = as.matrix(Xs)
  Ms = as.matrix(Ms)
  LA = length(unique(A))    
   if(LA<2) stop("Treatment levels are less than 2")

    RegGamma = glm(CC ~ Xs+A*Ms, family = "Gamma"(link="identity"),data = data);
    mus.gam = matrix(NA,N,LA)
    for(j in 1:LA){
      newdat = data.frame(Xs,A=rep(sort(unique(A))[j],N),Ms)
      mus.gam[,j] = predict(RegGamma,newdata=newdat)
    }

  out.gam = list(mus.gam, RegGamma);names(out.gam) = c("mus.gam","Gamma");out.gam
}


# Compute IPW estimates for individual treatment effect 
# Q: Observed effectiveness; C: observed cumulative cost outcome 
# A: Treatment; est_ps:estimated PS
# cen: Censoring indicator; est_cen: censoring probability 
# lambda: Willingness-to-pay parameter 
ITR_IPW = function(Q,C,A,est_ps,cen,est_cen,lambda){
  level.A = sort(unique(A))
  level.cen = sort(unique(cen))
  LA = length(level.A)
  if(LA<2) stop("Treatment levels are less than 2")
  if(dim(est_ps)[2]!=LA) stop("Different number of columns between treatment and PS")
  
  # IPW estimates
  muQ = muC = matrix(NA,N,LA)
  for(k in 1:LA){
    muQ[,k] = (cen/est_cen)*((A==level.A[k])*Q/est_ps[,k])
    muC[,k] = (cen/est_cen)*((A==level.A[k])*C/est_ps[,k])
  }

  # CE IPW estimates
  muY = lambda*muQ-muC
  
  # create contrast estimates
  Contrast.Y = Contrast.Q = Contrast.C = ITR_regCE = ITR_regST = rep(NA,N)
  for(j in 1:N){
    # CE 
    Contrast.Y[j] = max(muY[j,])-min(muY[j,])
    # Effectiveness (survival time)
    Contrast.Q[j] = max(muQ[j,])-min(muQ[j,])
    # Cumulative cost
    Contrast.C[j] = max(muC[j,])-min(muC[j,])
    
    # Treatment corresponds to preferred (higher) CE outcome
    ITR_regCE[j] = ifelse(cen[j]==0, NA, which(muY[j,]==max(muY[j,]))-1)
    # Treatment corresponds to preferred (higher) effectiveness outcome
    ITR_regST[j] = ifelse(cen[j]==0, NA, which(muQ[j,]==max(muQ[j,]))-1)
  }
  output = data.frame(Contrast.Q, Contrast.C, Contrast.Y, ITR_regCE, ITR_regST)
  output
}


# Compute IPW estimates for individual treatment effect 
# Q: Observed effectiveness; C: observed cumulative cost outcome 
# muQ.reg: Estimated effectiveness outcome; muC.reg: Estimated cumulative cost outcome
# A: Treatment; est_ps:estimated PS
# cen: Censoring indicator; est_cen: censoring probability 
# lambda: Willingness-to-pay parameter 
ITR_AIPW = function(Q,C,A,est_ps,cen,est_cen,muQ.reg,muC.reg,lambda){
  level.A = sort(unique(A))
  level.cen = sort(unique(cen))
  LA = length(level.A)
  N = length(A)
  if(LA<2) stop("Treatment levels are less than 2")
  if(dim(est_ps)[2]!=LA|dim(muQ.reg)[2]!=LA|dim(muC.reg)[2]!=LA) stop("Number of columns are inconsistent")
  
  # AIPW estimates
  muQ = muC = matrix(NA,N,LA)
  for(k in 1:LA){
    muQ[,k] = (cen/est_cen)*((A==level.A[k])*Q/est_ps[,k]+(1-(A==level.A[k])/est_ps[,k])*muQ.reg[,k])
    muC[,k] = (cen/est_cen)*((A==level.A[k])*C/est_ps[,k]+(1-(A==level.A[k])/est_ps[,k])*muC.reg[,k])
  }
  
  # CE AIPW estimates
  muY = lambda*muQ-muC
  
  # create contrast estimates
  Contrast.Y = Contrast.Q = Contrast.C = ITR_regCE = ITR_regST = rep(NA,N)
  for(j in 1:N){
    # CE 
    Contrast.Y[j] = max(muY[j,])-min(muY[j,])
    # Effectiveness (survival time)
    Contrast.Q[j] = max(muQ[j,])-min(muQ[j,])
    # Cumulative cost
    Contrast.C[j] = max(muC[j,])-min(muC[j,])
    
    # Treatment corresponds to preferred (higher) CE outcome
    ITR_regCE[j] = ifelse(cen[j]==0, NA, which(muY[j,]==max(muY[j,]))-1)
    # Treatment corresponds to preferred (higher) effectiveness outcome
    ITR_regST[j] = ifelse(cen[j]==0, NA, which(muQ[j,]==max(muQ[j,]))-1)
  }
  output = data.frame(Contrast.Q, Contrast.C, Contrast.Y, ITR_regCE, ITR_regST)
  output
}
