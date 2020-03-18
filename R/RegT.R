#---------------------------------------------------------------------------------------------------#######
#----------- Estimate survival time using Cox PH model                          --------------------#######
#----------- Outcome: ST, Follow-up time ST                                     --------------------#######
#----------- Event Indicator: event                                             --------------------#######
#----------- Predictor: treatment A, main covariates Xs, and interaction Ms(*A) --------------------#######
#---------------------------------------------------------------------------------------------------#######
library(survival);library(cubature)
Reg.mu.exp = function(ST,A,Xs,Ms,event,data){

  if(!(length(A) == dim(Xs)[1]))stop("Variables length do not match!")
  if(ncol(as.matrix(A))!=1 | length(unique(A))!= 2) stop("Multiple treatments present, or treatment levels are not 2!")
  if(!(unique(A)[1] %in% c(0,1)) | !(unique(A)[2] %in% c(0,1))) stop("Treatment levels are not 0 and 1")

    N = dim(data)[1]; Xs = as.matrix(Xs); Ms = as.matrix(Ms)
    fitsurvT = survreg(Surv(time = ST, event = event) ~ A*Ms+Xs, dist="lognormal",data=data)
    
    # estimate location parameter ui (on the log scale) for each individual & under each treament arm 
    log.survT = matrix(NA,N,2)
    for(j in 1:2){
      newdat = data.frame(Xs,A=rep(sort(unique(A))[j],N),Ms)
      log.survT[,j] = log(predict(fitsurvT,newdata=newdat,type = "response"))
    }
    
    mulog0 = log.survT[,1];mulog1 = log.survT[,2];sdlog=summary(fitsurvT)$scale
    mus.survRT = matrix(NA,N,2)
    for (i in 1:N){
      # survival function for each arm, same sd is used
      S_t0 = function(t){plnorm(t,meanlog=mulog0[i],sdlog=sdlog,lower.tail = FALSE)}
      S_t1 = function(t){plnorm(t,meanlog=mulog1[i],sdlog=sdlog,lower.tail = FALSE)}
      
      # integrate from 0 to tau
      mus.survRT[i,1] = cubintegrate(S_t0,0,tau,method = "cuhre",relTol=1e-03,absTol=1e-9)$integral
      mus.survRT[i,2] = cubintegrate(S_t1,0,tau,method = "cuhre",relTol=1e-03,absTol=1e-9)$integral
    };mus.survRT
}




