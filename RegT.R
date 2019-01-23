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
