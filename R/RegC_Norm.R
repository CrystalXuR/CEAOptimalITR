#---------------------------------------------------------------------------------------------------#######
#----------- Estimate cumulative cost using Lognormal model                     --------------------#######
#----------- Outcome: CC, cumulative cost                                       --------------------#######
#----------- Predictor: treatment A, main covariates Xs, and interaction Ms(*A) --------------------#######
#---------------------------------------------------------------------------------------------------#######

Reg.mu.norm = function(CC,A,Xs,Ms,data){
  
  if(!(length(A) == dim(Xs)[1]))stop("Variables length do not match!")
  if(ncol(as.matrix(A))!=1 | length(unique(A))!= 2) stop("Multiple treatments present, or treatment levels are not 2!")
  if(!(unique(A)[1] %in% c(0,1)) | !(unique(A)[2] %in% c(0,1))) stop("Treatment levels are not 0 and 1")
  
  RegNormal = glm(CC ~ A*Ms+Xs, family = "gaussian"(link="log"),data = data)
  newdat0 = data.frame(rep(0,dim(data)[1]),Ms,Xs);names(newdat0)[1] = "A"
  newdat1 = data.frame(rep(1,dim(data)[1]),Ms,Xs);names(newdat1)[1] = "A"
  mus.norm0 = predict(RegNormal,newdata=newdat0, type = "response")
  mus.norm1 = predict(RegNormal,newdata=newdat1, type = "response")
  mus.norm = data.frame(mus.norm0,mus.norm1);names(mus.norm) = c("muC0","muC1"); mus.norm 
}
