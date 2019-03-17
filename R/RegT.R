#---------------------------------------------------------------------------------------------------#######
#----------- Estimate survival time using Cox PH model                          --------------------#######
#----------- Outcome: ST, Follow-up time ST                                     --------------------#######
#----------- Event Indicator: event                                             --------------------#######
#----------- Predictor: treatment A, main covariates Xs, and interaction Ms(*A) --------------------#######
#---------------------------------------------------------------------------------------------------#######

Reg.mu.exp = function(ST,A,Xs,Ms,event,data){
  
  if(!(length(A) == dim(Xs)[1]))stop("Variables length do not match!")
  if(ncol(as.matrix(A))!=1 | length(unique(A))!= 2) stop("Multiple treatments present, or treatment levels are not 2!")
  if(!(unique(A)[1] %in% c(0,1)) | !(unique(A)[2] %in% c(0,1))) stop("Treatment levels are not 0 and 1")
    
    fitexpT = survreg(Surv(time = ST, event = event) ~ A*Ms+Xs, dist="exponential",data=data)
    newdat0 = data.frame(rep(0,dim(data)[1]),Ms,Xs);names(newdat0)[1] = "A"
    newdat1 = data.frame(rep(1,dim(data)[1]),Ms,Xs);names(newdat1)[1] = "A"
    survT0 = predict(fitexpT,newdata=newdat0, type = "response")
    survT1 = predict(fitexpT,newdata=newdat1, type = "response")
    mus.ExpT = data.frame(survT0,survT1);names(mus.ExpT) = c("muT0","muT1"); mus.ExpT     
}




  