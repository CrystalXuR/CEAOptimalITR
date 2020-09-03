#####--------------------------------------------------------------------------------#####
#####----- Functions Created for the Paper Entitled "Estimating the Optimal    ------#####
#####----- Individualized Treatment Rule from A Cost-Effectiveness Perspective ------#####
#####--------------------------------------------------------------------------------#####
library(survival);library(rpart);library(dplyr);library(cubature)
# Inputs: covariate matrix - X; covariates in main terms of outcome models- Xs;
#         covariates in interaction terms of outcome models- Ms; treatment - A; 
#         observed survival - ST; event indicator - event; observed total cost - C;
#         indicator of having event or being administratively censored - cen;
#         probability of uncensored (using cen as the response variable) - est_cen;
#         restriction time of the study - tau;willingness-to-pay threshold - lambda;
#         weight function - wf
OptimalCEITR <- function(X,Xs,Ms,A,ST,event,C,cen,est_cen,tau,lambda,wf="ITR_AIPW",data){
  
  # Function: Estimate Propensity Score 
  .EstPS <- function(A,X){
    
    if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) 
      stop("Treatment levels are not 0 and 1")
    
    dat <- data.frame(A,X)
    fit <- glm(A ~., family="binomial", data=dat)
    e   <- predict(fit,dat,type = "response")
    PS  <- data.frame(1-e,e);names(PS) = c("PrA0","PrA1")
    return(PS) 
  }
  
  # Function: Estimate Restricted Survival Time 
  .Reg.mu.surv <- function(ST,event,A,Xs,Ms,tau,data){
    
    N <- dim(data)[1]; Xs <- as.matrix(Xs); Ms <- as.matrix(Ms)
    
    if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) 
      stop("Treatment levels are not 0 and 1")
    
    fitsurvT <- survreg(Surv(time = ST, event = event) ~ Xs+A+A:Ms, dist="lognormal", data=data)
    
    # estimate location parameter ui (on the log scale) for each individual & under each treatment arm 
    log.survT <- matrix(NA,N,2)
    for(j in 1:2){
      newdat <- data; newdat$A <- rep(sort(unique(A))[j],N)
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
    
    outsurvRT <- list(mus.survRT, fitsurvT); names(outsurvRT) <- c("mus.survRT","SurvModel")
    return(mus.survRT)
  }
  
  
  # Function: Estimate cumulative cost 
  .Reg.mu.norm <- function(C,A,Xs,Ms,data){
    
    N <- dim(data)[1]; Xs <- as.matrix(Xs); Ms <- as.matrix(Ms)
    
    if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) 
      stop("Treatment levels are not 0 and 1")
    
    RegNormal <- glm(C ~ Xs+A+A:Ms,family = gaussian(link="log"),data = data)
    mus.norm  <- matrix(NA,N,2)
    for(j in 1:2){
      newdat <- data; newdat$A <- rep(sort(unique(A))[j],N)
      mus.norm[,j] <- predict(RegNormal,newdata=newdat,type="response")
    }
    
    out.norm <- list(mus.norm, RegNormal); names(out.norm) <- c("mus.norm","LogNormal")
    return(mus.norm)
  }
  
  
  # Function: Compute IPW estimates for classification weights 
  .ITR_IPW <- function(ST,C,A,est_ps,cen,est_cen,tau,lambda){
    
    N <- length(A)
    Q <- pmin(ST,rep(tau,N))
    if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) 
      stop("Treatment levels are not 0 and 1")
    
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
    
    # Optimal treatment w.r.t. NMB & effectiveness
    ITR_ipwCE <- ifelse(cen==0, NA, as.numeric(I(ipwY1-ipwY0>0)))
    
    output <- data.frame(delta.Y, ITR_ipwCE)
    return(output)
  }
  
  
  # Function: Compute AIPW estimates for classification weights 
  .ITR_AIPW <- function(Q,C,A,est_ps,cen,est_cen,muQ.reg,muC.reg,tau,lambda){
    
    N <- length(A)
    Q <- pmin(ST,rep(tau,N))
    if(sort(unique(A))[1]!=0 | sort(unique(A))[2]!=1) 
      stop("Treatment levels are not 0 and 1")
    
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
    
    # Optimal treatment w.r.t. NMB & effectiveness
    ITR_aipwCE <- ifelse(cen==0, NA, as.numeric(I((aipwY1-aipwY0)>0)))
    
    output <- data.frame(delta.Y, ITR_aipwCE)
    return(output)
  }
  
  PS <- .EstPS(A,X)
  RmuT <- .Reg.mu.surv(ST,event,A,Xs,Ms,tau,data)
  RmuC <- .Reg.mu.norm(C,A,Xs,Ms,data)
  if(wf=="ITR_IPW"){
    Out <- .ITR_IPW(Q,C,A,est_ps=PS,cen,est_cen,tau=tau,lambda=lambda)
  }else{
    Out <- .ITR_AIPW(Q,C,A,est_ps=PS,cen,est_cen,muQ.reg=RmuT,muC.reg=RmuC,tau=tau,lambda=lambda)
  }
  
  # Plug in classification weight (cw) and class label (cl) to decision trees
  cw <- Out$delta.Y
  cl <- Out$ITR_aipwCE
  tmpdat   <- data.frame(X,cl,cw)
  comp.dat <- tmpdat[complete.cases(tmpdat)==T,]
  cov <- comp.dat[,-c(dim(comp.dat)[2])]
  fit <- rpart(cl ~ ., weights=comp.dat$cw, method="class",data=cov)
  OptimalITR <- list(NA,2)
  OptimalITR[[1]] <- as.numeric(predict(fit, type="class"))-1
  OptimalITR[[2]] <- fit$variable.importance
  return(OptimalITR)
}


