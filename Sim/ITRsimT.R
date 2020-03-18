#####---------------------------------------------------------------------------#####
#####--------------------- Functions for simulating CE data --------------------#####
#####--------------------- Effect Modification on T only    --------------------#####

# Simulate survival time with exponential COX PH model 
surT.exp = function(i,trt,h0,betaF,gammaF,St,X0){
  tmpx = c(1,X0[i,1],X0[i,2],X0[i,3],X0[i,4],X0[i,5])                                                     
  ht = h0*exp(t(betaF)%*%tmpx-(gammaF[1]*X0[i,1]+gammaF[2]*X0[i,2])*trt[i])         # constant hazard rate 
  SurvT = -log(St[i])/ht                                                            # -log(S(t)) = cumulative hazard
  return(SurvT)
}

# simulate cumulative cost with gamma dist. 
cumC.gam = function(i, trt, k, betaC,gammaC,Nmonth,death,X0){
  tmpx = c(1,X0[i,1],X0[i,2],X0[i,3],X0[i,4],X0[i,5]) 
  tmp.beta = exp(t(betaC)%*%tmpx+gammaC*trt[i])
  
  inital.C  = rgamma(1, shape = k, scale = 5*tmp.beta)
  ongoing.C = rgamma(1, shape = k, scale = 0.5*tmp.beta)
  dying.C   = rgamma(1, shape = k, scale = 1.5*tmp.beta)
  
  cumC = inital.C+ongoing.C*Nmonth[i]+dying.C*death[i]
  return(cumC)
}
