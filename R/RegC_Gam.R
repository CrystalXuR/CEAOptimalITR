#---------------------------------------------------------------------------------------------------#######
#----------- Estimate cumulative cost using Gamma (identity) model              --------------------#######
#----------- Outcome: CC, cumulative cost                                       --------------------#######
#----------- Predictor: treatment A, main covariates Xs, and interaction Ms(*A) --------------------#######
#---------------------------------------------------------------------------------------------------#######

Reg.mu.gam = function(CC,A,Xs,Ms,event,data){
  
  if(!(length(A) == dim(Xs)[1]))stop("Variables length do not match!")
  if(ncol(as.matrix(A))!=1 | length(unique(A))!= 2) stop("Multiple treatments present, or treatment levels are not 2!")
  if(!(unique(A)[1] %in% c(0,1)) | !(unique(A)[2] %in% c(0,1))) stop("Treatment levels are not 0 and 1")
  
  RegGamma = glm(CC ~ A*Ms+Xs, family = "Gamma"(link="identity"), data = data)
  newdat0 = data.frame(rep(0,dim(data)[1]),Ms,Xs);names(newdat0)[1] = "A"
  newdat1 = data.frame(rep(1,dim(data)[1]),Ms,Xs);names(newdat1)[1] = "A"
  mus.gam0 = predict(RegGamma,newdata=newdat0, type = "response")
  mus.gam1 = predict(RegGamma,newdata=newdat1, type = "response")
  mus.gam = data.frame(mus.gam0,mus.gam1);names(mus.gam) = c("muC0","muC1"); mus.gam   
}
