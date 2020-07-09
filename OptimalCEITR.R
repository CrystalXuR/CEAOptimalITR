#####--------------------------------------------------------------------------------#####
#####----- Functions Created for the Paper Entitled "Estimating the Optimal    ------#####
#####----- Individualized Treatment Rule from A Cost-Effectiveness Perspective ------#####
#####--------------------------------------------------------------------------------#####

# Function: Estimate Propensity Score 
# Inputs: treatment - A; covariate matrix - X
EstPS <- function(A,X){
  
  if(length(A)!= nrow(as.matrix(X))) stop("Treatment and covariates have inconsistent dimensions")
  if(ncol(as.matrix(A))!=1 | length(unique(A))!=2) stop("Multiple treatments present or treatment is not binary")
  if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) stop("Treatment levels are not 0 and 1")
  
  dat <- data.frame(A,X)
  fit <- glm(A ~., family="binomial", data=dat)
  e   <- predict(fit,dat,type = "response")
  PS  <- data.frame(1-e,e);names(PS) = c("PrA0","PrA1");PS  
}

# Function: Estimate Restricted Survival Time 
# Inputs: outcome - ST; event indicator - event; treatment - A; covariates in main terms - Xs; covariates in interaction term - Ms
Reg.mu.surv <- function(ST,event,A,Xs,Ms,tau,data){
  
  N <- dim(data)[1]; Xs <- as.matrix(Xs); Ms <- as.matrix(Ms)
  
  if(length(ST)!=N | length(event)!=N | length(A)!=N | dim(Xs)[1]!=N | dim(Ms)[1]!=N) stop("Dimensions of input variables do not match")
  if(ncol(as.matrix(A))!=1 | length(unique(A))!=2) stop("Multiple treatments present or treatment is not binary")
  if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) stop("Treatment levels are not 0 and 1")
  
  fitsurvT <- survreg(Surv(time = ST, event = event) ~ Xs+A*Ms, dist="lognormal",data=data)
  
  # estimate location parameter ui (on the log scale) for each individual & under each treament arm 
  log.survT <- matrix(NA,N,2)
  for(j in 1:2){
    newdat <- data.frame(Xs,A=rep(sort(unique(A))[j],N),Ms)
    log.survT[,j] <- log(predict(fitsurvT,newdata=newdat,type = "response"))
  }
  
  mulog0 <- log.survT[,1]; mulog1 <- log.survT[,2]; sdlog <- summary(fitsurvT)$scale
  mus.survRT <- matrix(NA,N,2)
  for (i in 1:N){
    # survival function for each treatment arm
    S_t0 <- function(t){plnorm(t,meanlog=mulog0[i],sdlog=sdlog,lower.tail = FALSE)}
    S_t1 <- function(t){plnorm(t,meanlog=mulog1[i],sdlog=sdlog,lower.tail = FALSE)}
    # integrate from 0 to tau
    mus.survRT[i,1] <- cubintegrate(S_t0,0,tau,method = "cuhre",relTol=1e-04,absTol=1e-10)$integral
    mus.survRT[i,2] <- cubintegrate(S_t1,0,tau,method = "cuhre",relTol=1e-04,absTol=1e-10)$integral
  }
  
  outsurvRT <- list(mus.survRT, fitsurvT); names(outsurvRT) <- c("mus.survRT","SurvModel");outsurvRT
}

# Verifying the log(predict()) gives X*beta on log scale
# beta=summary(fitsurvT)$coefficients;BBeta=c(beta[1],beta[6:7],beta[2:5],beta[8:9])
# dat=as.matrix(data.frame(rep(1,N),data[,1:6],data$A*data$x1,data$A*data$x2))
# logmean = dat%*%BBeta
# mulog = log(predict(fitsurvT,type = "response")) # have verified this is beta%*%X, see above 
# RT=rep(NA,N)
# for (i in 1:N){
#   S_t = function(t){plnorm(t,meanlog=mulog[i],sdlog=sdlog,lower.tail = FALSE)}
#   RT[i]=integrate(S_t,0,tau)$value
# };summary(RT)


# Function: Estimate cumulative cost 
# Inputs: outcome - CC; treatment - A; covariates in main terms - Xs; covariates in interaction term: Ms
Reg.mu.norm <- function(CC,A,Xs,Ms,data){
  
  N <- dim(data)[1]; Xs <- as.matrix(Xs); Ms <- as.matrix(Ms)
  
  if(length(CC)!=N | length(A)!=N | dim(Xs)[1]!=N | dim(Ms)[1]!=N) stop("Dimensions of input variables do not match")
  if(ncol(as.matrix(A))!=1 | length(unique(A))!=2) stop("Multiple treatments present or treatment is not binary")
  if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) stop("Treatment levels are not 0 and 1")
  
  RegNormal <- glm(CC ~ Xs+A*as.matrix(Ms),family = "gaussian"(link="log"),data = data)
  mus.norm  <- matrix(NA,N,2)
  for(j in 1:2){
    newdat <- data.frame(Xs,A=rep(sort(unique(A))[j],N),Ms)
    mus.norm[,j] <- predict(RegNormal,newdata=newdat,type="response")
  }

  out.norm <- list(mus.norm, RegNormal); names(out.norm) <- c("mus.norm","LogNormal"); out.norm
}


# Function: Compute IPW estimates for classification weights 
# Inputs: observed effectiveness - Q; observed cost - C; treatment - A; estimated PS - est_ps; event indicator - cen;
#         probability of not being censored - est_cen; willingness-to-pay threshold - lambda 
ITR_IPW <- function(Q,C,A,est_ps,cen,est_cen,lambda){
  
  N <- length(A)
  if(length(Q)!=N | length(C)!=N | length(est_cen)!=N | length(cen)!=N | dim(est_ps)[1]!=N) stop("Dimensions of input variables do not match")
  if(ncol(as.matrix(A))!=1 | length(unique(A))!=2) stop("Multiple treatments present or treatment is not binary")
  if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) stop("Treatment levels are not 0 and 1")
  
  # IPW estimates for restricted survival time (Q) and cumulative cost (C)
  ipwQ0 <- ifelse(A==0, (cen/est_cen)*(Q/est_ps[,1]),0)
  ipwQ1 <- ifelse(A==1, (cen/est_cen)*(Q/est_ps[,2]),0)
  ipwC0 <- ifelse(A==0, (cen/est_cen)*(C/est_ps[,1]),0)
  ipwC1 <- ifelse(A==1, (cen/est_cen)*(C/est_ps[,2]),0)
  
  # IPW estimates for NMB
  ipwY0 <- lambda*ipwQ0-ipwC0
  ipwY1 <- lambda*ipwQ1-ipwC1
  
  # Estimated individual-level treatment effect 
  delta.Y <- abs(ipwY1-ipwY0)
  delta.Q <- abs(ipwQ1-ipwQ0)
  delta.C <- abs(ipwC1-ipwC0)
  
  # Optimal treatment w.r.t. NMB & effectiveness
  ITR_ipwCE <- ifelse(cen==0, NA, as.numeric(I(ipwY1-ipwY0>0)))
  ITR_ipwQ  <- ifelse(cen==0, NA, as.numeric(I(ipwQ1-ipwQ0>0)))
  
  output <- data.frame(delta.Y, delta.Q, delta.C, ITR_ipwCE, ITR_ipwQ);output
  
}


# Function: Compute AIPW estimates for classification weights 
# Inputs: observed effectiveness - Q; observed cost - C; treatment - A; estimated PS - est_ps; event indicator - cen;
#         probability of not being censored - est_cen; willingness-to-pay threshold - lambda 
#         estimated effectiveness outcome - muQ.reg; estimated cost outcome - muC.reg
ITR_AIPW <- function(Q,C,A,est_ps,cen,est_cen,muQ.reg,muC.reg,lambda){
  
  N <- length(A)
  if(length(Q)!=N | length(C)!=N | length(est_cen)!=N | length(cen)!=N | dim(est_ps)[1]!=N) stop("Dimensions of input variables do not match")
  if(ncol(as.matrix(A))!=1 | length(unique(A))!=2) stop("Multiple treatments present or treatment is not binary")
  if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) stop("Treatment levels are not 0 and 1")
  
  # AIPW estimates for restricted survival time (Q) and cumulative cost (C)
  aipwQ0 <- (cen/est_cen)*ifelse(A==0, Q/est_ps[,1]+(1-1/est_ps[,1])*muQ.reg[,1], muQ.reg[,1])
  aipwQ1 <- (cen/est_cen)*ifelse(A==1, Q/est_ps[,2]+(1-1/est_ps[,2])*muQ.reg[,2], muQ.reg[,2])
  aipwC0 <- (cen/est_cen)*ifelse(A==0, C/est_ps[,1]+(1-1/est_ps[,1])*muC.reg[,1], muC.reg[,1])
  aipwC1 <- (cen/est_cen)*ifelse(A==1, C/est_ps[,2]+(1-1/est_ps[,2])*muC.reg[,2], muC.reg[,2])
  
  # AIPW estimates for NMB
  aipwY0 <- lambda*aipwQ0-aipwC0
  aipwY1 <- lambda*aipwQ1-aipwC1
  
  # Estimated individual-level treatment effect 
  delta.Y <- abs(aipwY1-aipwY0)
  delta.Q <- abs(aipwQ1-aipwQ0)
  delta.C <- abs(aipwC1-aipwC0)
  
  # Optimal treatment w.r.t. NMB & effectiveness
  ITR_aipwCE <- ifelse(cen==0, NA, as.numeric(I((aipwY1-aipwY0)>0)))
  ITR_aipwQ  <- ifelse(cen==0, NA, as.numeric(I((aipwQ1-aipwQ0)>0)))
  
  output <- data.frame(delta.Y, delta.Q, delta.C, ITR_aipwCE, ITR_aipwQ);output
  
}
