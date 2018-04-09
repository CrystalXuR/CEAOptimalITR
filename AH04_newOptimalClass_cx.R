#sourceDirectory("C:/Users/u0918455/Box Sync/A_Crystal's Thesis/Rsource")

setGeneric(name = ".newOptimalClass_cx", 
           def = function(moPropen, moMain, moCensor, moCont, moClass, ...){
                   standardGeneric(".newOptimalClass_cx")
                 })

#----------------------------------------------------------------------#
# IPWE classification routine                                          #
#----------------------------------------------------------------------#
#  params                                                              #
#  moPropen : modelObj for propensity                                  #
#  moMain   : NULL                                                     #
#  moCont   : NULL                                                     #
#  moClass  : modelObj for classification                              #
#  data     : data.frame of covariates and treatment history           #
#  response : outcome of interest                                      #
#  txName   : character name of treatment in data                      #
#  iter     : ignored                                                  #
#  suppress : T/F indicating if screen prints are suppressed           #
#----------------------------------------------------------------------#
.ipwe <- function(moPropen, 
                  moMain,
                  moCensor,
                  moCont,
                  moClass, 
                  data, 
                  responseE,
                  responseC,
                  txName, 
                  iter,
                  suppress,
                  lambda,
                  censored,...) {

  if( !suppress ) {
    cat("Inverse Probability Weighted Estimator - Classification.\n")
  }

  #------------------------------------------------------------------#
  # Process treatment information.                                   #
  #------------------------------------------------------------------#
  txInfo <- .newTxInfo(fSet = NULL, 
                       txName = txName, 
                       data = data,
                       suppress = suppress)

  #------------------------------------------------------------------#
  # Fit propensity models                                            #
  #------------------------------------------------------------------#
  propen <- .newPropensityRegression(moPropen = moPropen, 
                                     txInfo = txInfo, 
                                     data = data,
                                     suppress = suppress)
  #------------------------------------------------------------------#
  # Fit ceonsoring weight models                                     #
  #------------------------------------------------------------------#
  Cweight <- .newPropensityRegression(moCensor = moCensor, 
                                      txInfo = txInfo, 
                                      data = data,
                                      suppress = suppress)
  #------------------------------------------------------------------#
  # Calculate contrast                                               #
  #------------------------------------------------------------------#
  contrast <- .contrastIPWE(txInfo = txInfo, 
                            propensity = propen, 
                            data = data,
                            Cweight = Cweight,
                            lambda = lambda,
                            censored = censored,
                            responseE = responseE,
                            responseC = responseC)

  #------------------------------------------------------------------#
  # Obtain classification fit                                        #
  #------------------------------------------------------------------#
  classFit <- .classFun(contrast = contrast$contrast, 
                        moClass = moClass,  
                        data = data)

  #------------------------------------------------------------------#
  # Obtain classification prediction from fit                        #
  #------------------------------------------------------------------#
  grp1 <- classFit$opt == "1"
  estResponse <- sum(contrast$contrast[grp1]) / nrow(data) + 
                 contrast$mean.mu0

  #------------------------------------------------------------------#
  # Estimate optimal treatment                                       #
  #------------------------------------------------------------------#
  optTx <- drop(predict(classFit$cf, newdata = data))

  if( !suppress ) {
    cat("Classification Analysis\n")
    print(classFit$cf)
    cat("Value:", estResponse, "\n")
  }

  oc1 <- new("OptimalClassIPWE",
             "classif"        = classFit$cf,
             "propen"         = propen,
             "estimatedValue" = estResponse,
             "optimalTx"      = optTx,
             "call"           = NULL)

  return( oc1 )

}

setMethod(f = ".newOptimalClass_cx",    
          signature = c(moPropen = "modelObj",
                        moMain   = "NULL",
                        moCensor = "modelObj",
                        moCont   = "NULL",
                        moClass  = "modelObj"), 
          definition = .ipwe)

#----------------------------------------------------------------------#
# prepares data and calls classification method specified by user.     #
#----------------------------------------------------------------------#
# contrast : Vector of the contrast function for each sample.          #
# moClass  : an object of class modelObj, which defines the models and #
#            R methods to be used to obtain parameter estimates and    #
#            predictions for the classification                        #
#            See ?modelObj for details.                                #
#                                                                      #
#            It is assumed that the solver.method contains             #
#              weights : A vector of weights to be used in the fitting #
#                        process.                                      #
# data     : data frame of covariates and response                     #
#----------------------------------------------------------------------#
#= Returns a list                                                     =#
#=   cf: classification fit object                                    =#
#=  opt: optimal treatment regime for training set                    =#
#----------------------------------------------------------------------#
.classFun <- function(contrast, moClass, data){

  #------------------------------------------------------------------#
  # Classification weight variable.                                  #
  #------------------------------------------------------------------#
  weights <- abs(contrast)

  #------------------------------------------------------------------#
  # Normalize weights                                                #
  #------------------------------------------------------------------#
  norm.weights <- weights / sum(weights)

  #------------------------------------------------------------------#
  # Add weights to formal arguments of classification method         #
  #------------------------------------------------------------------#
  cArgs <- solverArgs(moClass)
  cArgs[[ "weights" ]] <- norm.weights
  solverArgs(moClass) <- cArgs

  #------------------------------------------------------------------#
  # Classification labels                                            #
  #------------------------------------------------------------------#
  ZinternalZ <- as.numeric(contrast > -1.5e-8)
  ZinternalZ <- as.factor(ZinternalZ)

  #------------------------------------------------------------------#
  # Obtain classification fit using fit routine of modelObj          #
  #------------------------------------------------------------------#
  cf <- fit(object = moClass, data = data, response = ZinternalZ)

  #------------------------------------------------------------------#
  # Predict classification of training set                           #
  #------------------------------------------------------------------#
  pred <- predict(object = cf, newdata = data)

  if( !is(pred, "factor") ) {
    pred <- as.factor(pred)
  }
  levs <- levels(pred)

  if( length(levs) == 1L ) {
    tst <- c("0","1") %in% levs
    if( !any(tst) )  {
      msg <- "unable to resolve class designation returned by predict"
      UserError("input", msg)
    }
  } else if( length(levs) == 2L ) {
    tst <- c("0","1") %in% levs
    if( !all(tst) ) {
      pred <- factor(pred, labels=c("0","1"))
    }
  } else {
    msg <- "classification prediction method returns more than 2 classes"
    UserError("input", msg)
  }
  
  return(list("cf" = cf, "opt" = pred))
}


#----------------------------------------------------------------------#
# IPWE contrast function for a single decision point binary tx.        #
#----------------------------------------------------------------------#
# txInfo      : an object of class txInfo                              #
# propensity  : propensity fit object                                  #
# Cweight.trt : censoring weight in the treated group                  #
# Cweight.ctrl: censoring weight in the control group                  #
# data        : data frame of covariates                               #
# response    : a response vector                                      #
#----------------------------------------------------------------------#
#= Returns a list                                                     =#
#=    constrast, mean.mu0                                             =#
#----------------------------------------------------------------------#
.contrastIPWE <- function(txInfo, 
                          propensity,
                          Cweight,
                          data,
                          lambda,
                          censored,
                          responseE,
                          responseC){

  #------------------------------------------------------------------#
  # Obtain matrix of propensity for treatment                        #
  #------------------------------------------------------------------#
  pm <- predict(object = propensity, newdata = data)
  
  #------------------------------------------------------------------#
  # Retrieve treatment column header                                 #
  #------------------------------------------------------------------#
  txName <- .getTxName(txInfo)

  #------------------------------------------------------------------#
  # Retrieve super set of treatment options                          #
  #------------------------------------------------------------------#
  fSetSuperSet <- .getSuperSet(txInfo)

  #------------------------------------------------------------------#
  # Convert treatment to 0/1 notation                                #
  #------------------------------------------------------------------#
  txVec <- numeric(nrow(data))
  txVec[data[,txName] == fSetSuperSet[1L]] <- 0L
  txVec[data[,txName] == fSetSuperSet[2L]] <- 1L

  #------------------------------------------------------------------#
  # Obtain matrix of ceonsoring for each treatment group             #
  #------------------------------------------------------------------#
  trtdata = data
  trtdata[,txName]=rep(1,length(trtdata[,txName]))
  ctrldata = data
  ctrldata[,txName]=rep(0,length(ctrldata[,txName]))
  cw.trt <- predict(object = Cweight, newdata = trtdata)
  cw.ctrl <- predict(object = Cweight, newdata = ctrldata)
  
  #------------------------------------------------------------------#
  # Calculate IPWE contrast function.                                #
  #------------------------------------------------------------------#
  ymE <- { txVec / pm[,2L] - (1.0 - txVec) / pm[,1L] } * responseE
  ymC <- { txVec*censored / (pm[,2L]*cw.trt) - (1.0 - txVec)*censored / (pm[,1L]*cw.ctrl) } * responseC
  ym <- lambda*ymE-ymC

  #------------------------------------------------------------------#
  # Calculate the mean value of patients that received base tx       #
  #------------------------------------------------------------------#
  mmu <- lambda*(1.0 - txVec) / pm[,1L] * responseE - (1.0 - txVec) / (pm[,1L]*cw.ctrl) * responseC 
  mmu <- mean(mmu)

  return(list("contrast" = ym,
              "mean.mu0" = mmu))
}



