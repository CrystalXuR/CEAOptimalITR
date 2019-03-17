#----------------------------------------------------------------------------------------#######
#----------- Estimating treatment effect & treatment regime using IPW      --------------#######
#----------- Input: follow-up time (ST), cumulative cost (CC),             --------------#######
#-----------        treatment (A), propensity score (est_ps),              --------------#######
#-----------        censoring indicator (cen), censoring weight (est_cen), --------------#######
#-----------        Willingness-to-pay parameter (lambda)                  --------------#######
#----------------------------------------------------------------------------------------#######

ITR_IPW = function(ST,CC,A,est_ps,cen,est_cen,lambda){
  
  if(!(length(A) == length(cen)))stop("Variables length do not match!")
  if(ncol(as.matrix(A))!=1 | length(unique(A))!= 2) stop("Multiple treatments present, or treatment levels are not 2!")
  if(!(unique(A)[1] %in% c(0,1)) | !(unique(A)[2] %in% c(0,1))) stop("Treatment levels are not 0 and 1")
  
  # IPW estimates for survival time (ST) and cumulative cost
  ipwST0 = ifelse(A==0, (cen/est_cen)*(ST/est_ps[,1]),0)
  ipwST1 = ifelse(A==1, (cen/est_cen)*(ST/est_ps[,2]),0)
  ipwCC0 = ifelse(A==0, (cen/est_cen)*(CC/est_ps[,1]),0)
  ipwCC1 = ifelse(A==1, (cen/est_cen)*(CC/est_ps[,2]),0)
  
  # CE IPW estimates
  ipwY0 = lambda*ipwST0-ipwCC0
  ipwY1 = lambda*ipwST1-ipwCC1
  
  # Estimated treatment effect & treatment regime
  for(j in 1:N){

    delta.Y = abs(ipwY1-ipwY0)
    delta.ST = abs(ipwST1-ipwST0)
    delta.CC = abs(ipwCC1-ipwCC0)
    
    # Treatment corresponds to preferred (higher value) outcomes
    ITR_ipwCE = ifelse(cen==0, NA, as.numeric(I((ipwY1-ipwY0)>0)))
    ITR_ipwST = ifelse(cen==0, NA, as.numeric(I((ipwST1-ipwST0)>0)))
  }
  output = data.frame(delta.Y, delta.ST, delta.CC, ITR_ipwCE, ITR_ipwST);output
}



