#----------------------------------------------------------------------------------------#######
#----------- Estimating treatment effect & treatment regime using AIPW     --------------#######
#----------- Input: follow-up time (ST), cumulative cost (CC),             --------------#######
#-----------        estimated survival time (muT.reg),                     --------------#######
#-----------        estimated cumulative cost (muC.reg),                   --------------#######
#-----------        treatment (A), propensity score (est_ps),              --------------#######
#-----------        censoring indicator (cen), censoring weight (est_cen), --------------#######
#-----------        Willingness-to-pay parameter (lambda)                  --------------#######
#----------------------------------------------------------------------------------------#######

ITR_AIPW = function(ST,CC,A,est_ps,cen,est_cen,muT.reg,muC.reg,lambda){
  
  if(!(length(A) == length(cen)))stop("Variables length do not match!")
  if(ncol(as.matrix(A))!=1 | length(unique(A))!= 2) stop("Multiple treatments present, or treatment levels are not 2!")
  if(!(unique(A)[1] %in% c(0,1)) | !(unique(A)[2] %in% c(0,1))) stop("Treatment levels are not 0 and 1")
  
  # IPW estimates for survival time (ST) and cumulative cost
  aipwST0 = (cen/est_cen)*ifelse(A==0, ST/est_ps[,1], muT.reg[,1]/est_ps[,1])
  aipwST1 = (cen/est_cen)*ifelse(A==1, ST/est_ps[,2], muT.reg[,2]/est_ps[,2])
  aipwCC0 = (cen/est_cen)*ifelse(A==0, CC/est_ps[,1], muC.reg[,1]/est_ps[,1])
  aipwCC1 = (cen/est_cen)*ifelse(A==1, CC/est_ps[,2], muC.reg[,2]/est_ps[,2])
  
  # CE IPW estimates
  aipwY0 = lambda*aipwST0-aipwCC0
  aipwY1 = lambda*aipwST1-aipwCC1

  # Estimated treatment effect & treatment regime
  for(j in 1:N){
    
    delta.Y = abs(aipwY1-aipwY0)
    delta.ST = abs(aipwST1-aipwST0)
    delta.CC = abs(aipwCC1-aipwCC0)
    
    # Treatment corresponds to preferred (higher value) outcomes
    ITR_aipwCE = ifelse(cen==0, NA, as.numeric(I((aipwY1-aipwY0)>0)))
    ITR_aipwST = ifelse(cen==0, NA, as.numeric(I((aipwST1-aipwST0)>0)))
  }
  output = data.frame(delta.Y, delta.ST, delta.CC, ITR_aipwCE, ITR_aipwST);output
}