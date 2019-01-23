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
