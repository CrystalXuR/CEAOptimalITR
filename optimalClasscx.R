#### Optimal Class Source Code ####
source("C:/Users/u0918455/Box Sync/A_Crystal's Thesis/1 - Optimization-Classification/2 - Simulations & implementation/Simulation Codes and Results/Optimal Class Code/AH04_newOptimalClass_cx.r")

.optimalClasscx<-function (..., moPropen, moMain, moCensor, moCont, moClass, 
          lambda, censored, data, responseE, responseC,
          txName, iter = 0L, verbose = TRUE) 
{
  if (missing(moPropen)) 
    moPropen <- NULL
  if (is(moPropen, "NULL")) {
    UserError("input", "moPropen must be provided")
  }
  else if (!is(moPropen, "modelObj")) {
    UserError("input", "moPropen must be an object of class modelObj")
  }
  if (missing(moMain)) 
    moMain <- NULL
  if (!is(moMain, "modelObj") && !is(moMain, "NULL")) {
    UserError("input", "moMain must be an object of class modelObj or NULL")
  }
  if (missing(moCensor)) 
    moCensor <- NULL
  if (is(moCensor, "NULL")) {
    UserError("input", "moCensor must be provided")
  }
  else if (!is(moCensor, "modelObj")) {
    UserError("input", "moCensor must be an object of class modelObj")
  }
  if (missing(moCont)) 
    moCont <- NULL
  if (!is(moCont, "modelObj") && !is(moCont, "NULL")) {
    UserError("input", "moCont must be an object of class modelObj or NULL")
  }
  if (missing(moClass)) 
    moClass <- NULL
  if (is(moClass, "NULL")) {
    UserError("input", "moClass must be provided")
  }
  if (!is(moClass, "modelObj")) {
    UserError("input", "moClass must be an object of class modelObj")
  }
  if (!is(data, "data.frame")) {
    UserError("input", "'data' must be a data.frame")
  }
  if (is(responseE, "data.frame")) 
    responseE <- data.matrix(responseE)
  if (!is(responseE, "matrix")) {
    responseE <- matrix(data = responseE, ncol = 1L)
  }
  if (ncol(responseE) != 1L) {
    UserError("input", "'responseE' must be a vector")
  }
  if (is(responseC, "data.frame")) 
    responseC <- data.matrix(responseC)
  if (!is(responseC, "matrix")) {
    responseC <- matrix(data = responseC, ncol = 1L)
  }
  if (ncol(responseC) != 1L) {
    UserError("input", "'responseC' must be a vector")
  }
  responseE <- drop(responseE)
  responseC <- drop(responseC)
  data <- .checkTxData(txName, data)
  txVec <- .checkBinaryTx(txName, data)
  tx <- txVec
  tx[txVec < 0] <- 0L
  tx[txVec > 0] <- 1L
  if (!isTRUE(all.equal(tx, data[, txName]))) {
    cat("Treatment variable converted to {0,1}\n")
    data[, txName] <- as.integer(round(tx, 0))
  }
  if (!is(iter, "integer")) 
    iter <- as.integer(round(iter, 0L))
  if (iter < 0) 
    iter <- 0L
  if (!is(verbose, "logical")) {
    UserError("input", "'verbose' must be a TRUE/FALSE")
  }
  result <- .newOptimalClass_cx(moPropen = moPropen, moMain = moMain, moCensor = moCensor,lambda = lambda, censored = censored,
                                moCont = moCont, moClass = moClass, data = data, responseE = responseE, responseC = responseC, 
                                txName = txName, iter = iter, suppress = !verbose)
  result@call <- match.call()
  return(result)
}