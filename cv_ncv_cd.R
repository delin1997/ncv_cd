source("ncv_cd.R")

cv_ncv_cd <- function(X, y, lambda.min, lambda, nlambda=100, penalty=c("MCP", "SCAD", "lasso"), gamma=switch(penalty, SCAD=3.7, MCP=3), nfolds=10, eps=1e-4){
  X <- scale(X)
  y <- scale(y, scale = FALSE)
  n <- nrow(X)
  p <- ncol(X)
  if(missing(lambda)){
    if(missing(lambda.min)){
      lambda.min <- ifelse(n>p, .001, .05)
    }
    lambda.max <- max(t(X)%*%y/n)
    if(lambda.min==0){
      lambda <- exp(seq(log(lambda.min), log(0.001*lambda.max), length.out = nlambda))
    }else{
      lambda <- exp(seq(log(lambda.min), log(lambda.min*lambda.max), length.out = nlambda))
    }
  }
  MSE <- matrix(0, nrow = length(lambda), ncol = nfolds)
  index <- sample(1:nfolds, replace = T, size = n)
  for (i in 1:nfolds) {
    mod <- ncv_cd(X = X[index!=i, ], y = y[index!=i], lambda = lambda, penalty = penalty, gamma = gamma, eps = eps)
    beta <- mod$beta
    MSE[, i] <- apply(y[index==i] - X[index==i, ]%*%beta, 2, norm, type = "2")
  }
  w <- apply(matrix(1:nfolds), 1, function(x){
    sum(index==x)/n
  })
  MSE_aver <- MSE%*%w
  lambda_selected <- lambda[which.min(MSE_aver)]
  mod <- ncv_cd(X = X, y = y, lambda = lambda, penalty = penalty, gamma = gamma, eps = eps)
  return(list(beta = mod$beta, lambda = lambda, beta_cv = mod$beta[, which.min(MSE_aver)], lambda_cv = lambda_selected))
}