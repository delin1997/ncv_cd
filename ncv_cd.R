S <- function(z, lambda){
  sign(z)*max(0, abs(z)-lambda)
}


MCP <- function(z, lambda, gamma){
  if(abs(z)<=gamma*lambda){
    S(z,lambda)/(1-1/gamma)
  }else{
    z
  }
}

SCAD <- function(z, lambda, gamma){
  if(abs(z)<=2*lambda){
    S(z, lambda)
  }else if(abs(z)<=gamma*lambda){
    S(z, gamma*lambda/(gamma-1))/(1-1/(gamma-1))
  }else{
    z
  }
}

ncv_cd <- function(X, y, family=c("gaussian", "binomial"), lambda.min, lambda, nlambda=100, penalty=c("MCP", "SCAD", "lasso"), gamma=switch(penalty, SCAD=3.7, MCP=3)){
  X <- scale(X)
  y <- scale(y, scale = FALSE)
  n <- nrow(X)
  p <- ncol(X)
  epsilon <- 1e-2*sqrt(p)
  if(missing(lambda)){
    if(missing(lambda.min)){
      lambda.min <- ifelse(n>p, .001, .05)
    }
    lambda.max <- max(t(X)%*%y/n)
    lambda <- exp(seq(log(lambda.min), log(lambda.max), length.out = nlambda))
  }
  beta_lambda <- matrix(0, nrow = p, ncol = length(lambda))
  beta <- rep(0, p) -> beta_new
  z <- rep(0, p)
  for (i in 1:length(lambda)) {
    r <- y - X%*%beta
    t <- 0
    repeat{
      for (j in 1:p) {
        z[j] <- X[, j]%*%r/n + beta[j]
        beta_new[j] <- get(penalty)(z = z[j], lambda = lambda[i], gamma = gamma)
        r <- r - (beta_new[j]-beta[j])*X[, j]
      }
      t <- t+1
      if(norm(beta_new-beta,"2")<epsilon){break}
      beta <- beta_new
    }
    beta_lambda[, i] <- beta
  }
  return(list(beta = beta_lambda, lambda = lambda))
}