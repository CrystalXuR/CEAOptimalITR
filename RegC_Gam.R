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