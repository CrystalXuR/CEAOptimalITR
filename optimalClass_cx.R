source("C:/Users/u0918455/Box Sync/A_Crystal's Thesis/1 - Optimization-Classification/2 - Simulations & implementation/Simulation Codes and Results/Optimal Class Code/AH04_newOptimalClass_cx.r")

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                      #
# implementation of estimation of optimal DTRs from a classification   #
# perspective.                                                         #
#                                                                      #
#  Baqun Zhang, Anastasios A. Tsiatis, Marie Davidian, Min Zhang and   #
#  Eric B. Laber. "Estimating optimal tx regimes from a classification #
#  perspective." Stat 2012; 1: 103-114.                                #
#                                                                      #
# Note that this method is a single decision point, binary treatment   #
# method                                                               #
#                                                                      #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# moPropen: an object of class modelObj, which defines the models and R#
#           methods to be used to obtain parameter estimates and       #
#           predictions for the propensity for treatment.              #
#           See ?modelObj for details.                                 #
#                                                                      #
#           If the prediction method specified in moPropen returns     #
#           predictions for only a subset of the categorical tx data,  #
#           it is assumed that the base level defined by levels(tx)[1] #
#           is the missing category.                                   #
#                                                                      #
# moMain  : an object of class modelObj, which defines the models and R#
#           methods to be used to obtain parameter estimates and       #
#           predictions for for the main effects component of the      #
#           outcome regression.                                        #
#           See ?modelObj for details.                                 #
#           NULL is an appropriate value.                              #
#                                                                      #
# moCont  : an object of class modelObj, which defines the models and R#
#           methods to be used to obtain parameter estimates and       #
#           predictions for for the contrasts component of the         #
#           outcome regression.                                        #
#           See ?modelObj for details.                                 #
#           NULL is an appropriate value.                              #
#                                                                      #
# moClass : an object of class modelObj, which defines the             #
#           models and R methods to be used to obtain parameter        #
#           estimates and predictions for the classification.          #
#           See ?modelObj for details.                                 #
#                                                                      #
# data    : a data frame of the covariates and tx histories            #
#           tx variable will be recast as factor if not provided as    #
#           such.                                                      #
#                                                                      #
# response: response vector                                            #
#                                                                      #
# txName  : an character giving the column header of the column in data#
#           that contains the tx covariate.                            #
#                                                                      #
# iter    : an integer                                                 #
#                                                                      #
#           >=1 if moMain and moCont are to be fitted iteratively      #
#           The value is the maximum number of iterations.             #
#           Note the iterative algorithms is as follows:               #
#           Y = Ymain + Ycont                                          #
#            (1) hat(Ycont) = 0                                        #
#            (2) Ymain = Y - hat(Ycont)                                #
#            (3) fit Ymain ~ moMain                                    #
#            (4) set Ycont = Y - hat(Ymain)                            #
#            (5) fit Ycont ~ moCont                                    #
#            (6) Repeat steps (2) - (5) until convergence or           #
#            a maximum of iter iterations.                             #
#                                                                      #
#           <=0 moMain and moCont will be combined and fit as a single #
#           object.                                                    #
#                                                                      #
#           Either categorical or integer data can be provided for the #
#           tx. If categorical, the fitted contrast and main effects   #
#           are defined relative to the base category {defined as      #
#           levels()[1]}. The values may not be those returned by      #
#           predict(object) if iterate fits are used. If integer, the  #
#           fitted contrast and main effects are defined               #
#           relative to no tx (tx = 0).                                #
#                                                                      #
#           Note that if iter <= 0, all non-model components of the    #
#           moMain and moCont must be identical                        #
#                                                                      #
# verbose:  a logical                                                  #
#           If FALSE, screen prints are suppressed.                    #
#                                                                      #
#======================================================================#
#=                                                                    =#
#= Returns an object of class optimalClass                            =#
#=                                                                    =#
#======================================================================#
optimalClasscx <- function(..., 
                         moPropen,
                         moMain,
                         moCensor,
                         moCont,
                         moClass,
                         lambda,
                         censored,
                         data,
                         responseE,
                         responseC,
                         txName,
                         iter = 0L,
                         verbose = TRUE){

  #------------------------------------------------------------------#
  # Ensure that moPropen is modelObj                                 #
  #------------------------------------------------------------------#
  if( missing(moPropen) ) moPropen <- NULL
  if( is(moPropen, "NULL") ) {
    UserError("input", "moPropen must be provided")
  } else if( !is(moPropen, "modelObj") ) {
    UserError("input",
              "moPropen must be an object of class modelObj")
  }

  #------------------------------------------------------------------#
  # Ensure that moMain is modelObj or NULL                           #
  #------------------------------------------------------------------#
  if( missing(moMain) ) moMain <- NULL
  if( !is(moMain, "modelObj") && !is(moMain, "NULL") ) {
    UserError("input",
              "moMain must be an object of class modelObj or NULL")
  }

  #------------------------------------------------------------------#
  # Ensure that moCensor is modelObj                                 #
  #------------------------------------------------------------------#
  if( missing(moCensor) ) moCensor <- NULL
  if( is(moCensor, "NULL") ) {
    UserError("input", "moCensor must be provided")
  } else if( !is(moCensor, "modelObj") ) {
    UserError("input",
              "moCensor must be an object of class modelObj")
  }
  
  #------------------------------------------------------------------#
  # Ensure that moCont is modelObj or NULL.                          #
  #------------------------------------------------------------------#
  if( missing(moCont) ) moCont <- NULL
  if( !is(moCont, "modelObj") && !is(moCont, "NULL") ) {
    UserError("input",
              "moCont must be an object of class modelObj or NULL")
  }

  #------------------------------------------------------------------#
  # Ensure that a classification modelObj is provided.               #
  #------------------------------------------------------------------#
  if( missing(moClass) ) moClass <- NULL
  if( is(moClass, "NULL") ) {
    UserError("input", "moClass must be provided")
  }
  if( !is(moClass, "modelObj") ) {
    UserError("input",
              "moClass must be an object of class modelObj")
  }

  #------------------------------------------------------------------#
  # data must be provided as a data.frame object.                    #
  #------------------------------------------------------------------#
  if( !is(data, "data.frame") ) {
    UserError("input", "'data' must be a data.frame")
  }

  #------------------------------------------------------------------#
  # responseE must be vector                                          #
  #------------------------------------------------------------------#
  if( is(responseE, "data.frame") ) responseE <- data.matrix(responseE)
  if( !is(responseE, "matrix") ) {
    responseE <- matrix(data = responseE, ncol = 1L)
  }

  if( ncol(responseE) != 1L ) {
    UserError("input",
              "'responseE' must be a vector")
  }

  responseE <- drop(responseE)

  #------------------------------------------------------------------#
  # responseC must be vector                                          #
  #------------------------------------------------------------------#
  if( is(responseC, "data.frame") ) responseC <- data.matrix(responseC)
  if( !is(responseC, "matrix") ) {
    responseC <- matrix(data = responseC, ncol = 1L)
  }
  
  if( ncol(responseC) != 1L ) {
    UserError("input",
              "'responseC' must be a vector")
  }
  
  responseC <- drop(responseC)
  
  
  #------------------------------------------------------------------#
  # Verify treatment is appropriately coded.                         #
  #------------------------------------------------------------------#
  data <- .checkTxData(txName, data)

  #------------------------------------------------------------------#
  # Treatments must be binary                                        #
  # Note that NAs are allowed                                        #
  #------------------------------------------------------------------#
  txVec <- .checkBinaryTx(txName, data)
  tx <- txVec
  tx[txVec < 0] <- 0L
  tx[txVec > 0] <- 1L

  if( !isTRUE(all.equal(tx, data[,txName])) ) {
    cat("Treatment variable converted to {0,1}\n")
    data[,txName] <- as.integer(round(tx,0))
  }

  if( !is(iter, "integer") ) iter <- as.integer(round(iter,0L))
  if( iter < 0 ) iter <- 0L

  if( !is(verbose, "logical") ) {
    UserError("input",
              "'verbose' must be a TRUE/FALSE")
  }

  result <- .newOptimalClass_cx(moPropen = moPropen,
                                moMain = moMain,
                                moCensor = moCensor,
                                moCont = moCont,
                                moClass = moClass,
                                lambda = lambda, 
                                censored = censored,
                                data = data,
                                responseE = responseE,
                                responseC = responseC,
                                txName = txName,
                                iter = iter,
                                suppress = !verbose)

  result@call <- match.call()

  return(result)


}
