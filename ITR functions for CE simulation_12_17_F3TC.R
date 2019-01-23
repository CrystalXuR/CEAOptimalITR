#####---------------------------------------------------------------------------#####
#####--------------------- Functions for simulating CE data --------------------#####
#####------------------- Effect Modification on both T and C -------------------#####

# Simulate survival time with exponential COX PH model  
surT.exp = function(i,trt_i,h0,betaF,gammaF){
  tmpx = c(1, x1[i],x2[i],x3[i],x4[i],x5[i])                                                     
  haz.rate = h0 * exp(t(betaF) %*% tmpx - (gammaF[2]*x1[i]+gammaF[3]*x2[i])*trt_i[i])       # hazard at time t, does not depend on time t since h0 is constant
  return(-log(surv.prob[i])/(haz.rate))                                                     # -log(S(t)) = cumulative hazard. Survival T = Cumulative hazard/hazard at time t
}

# cumulative summation tool
sumfun = function(x,start,end){
  return(sum(x[start:end]))
}

# simulate cumulative cost with gamma dist. 
cumC.gam = function(i, trt_i, k, betaC,gammaC,I.delta){
  tmpx = c(1,x1[i],x2[i],x3[i],x4[i],x5[i]) 
  tmp.beta = exp(t(betaC) %*% tmpx + (gammaC[2]*x1[i]+gammaC[3]*x2[i])*trt_i[i])
  
  inital.C = rgamma(1,shape = k,scale=5*tmp.beta)
  ongoing.C = rgamma(48, shape = k, scale=0.5*tmp.beta)
  ongoing.C = I.delta[i,]*ongoing.C
  ongoing.C = ifelse(ongoing.C==0, NA,ongoing.C)
  dying.C = rgamma(1, shape = k, scale=1.5*tmp.beta)                    # if one is censored, we don't add the dying cost(=0)
  
  cumC = rep(NA,48)
  for (j in 1:48){
    cost.series = c(inital.C,ongoing.C)
    cumC[j] = sumfun(cost.series, 1,j)
  }
  
  ni = sum(!is.na(cumC))
  cumC[ni] = cumC[ni]+dying.C
  cumCost = cumC
  return(cumCost)
}

# use only one legend for plots
g_legend = function(a.gplot){
  tmp  =  ggplot_gtable(ggplot_build(a.gplot))
  leg  =  which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend  =  tmp$grobs[[leg]]
  return(legend)}

