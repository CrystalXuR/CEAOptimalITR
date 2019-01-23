# Estimate PS using logistic regression
# Outcome: treatment A 
# Predictor: covariate Xs
EstPS = function(A,Xs){
  if(ncol(as.matrix(A))!=1) stop("Multiple treatments are not allowed")
  if(length(A)!= nrow(as.matrix(Xs))) stop("Non-conformable between outcome and predictors")
  if(length(unique(A))<=1) stop("Outcome levels are less than 2")
  level.A = sort(unique(A))
  
  dat = data.frame(A,Xs)
  PSfit = glm(A ~., family="binomial", data=dat)
  p = predict(PSfit,dat,type = "response")
  PS = data.frame((1-p),p);names(PS) = paste("A=",level.A,sep="");PS  
}
