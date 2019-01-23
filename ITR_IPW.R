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