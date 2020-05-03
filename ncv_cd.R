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

ncv_cd <- function(X, y, family=c("gaussian", "binomial"), lambda.min, lambda, nlambda=100, penalty=c("MCP", "SCAD", "lasso"), gamma=switch(penalty, SCAD=3.7, MCP=3), plot = FALSE){
  X <- scale(X)
  y <- scale(y, scale = FALSE)
  n <- nrow(X)
  p <- ncol(X)
  epsilon <- 1e-4*sqrt(p)
  if(missing(lambda)){
    if(missing(lambda.min)){
      lambda.min <- ifelse(n>p, .001, .05)
    }
    lambda.max <- max(t(X)%*%y/n)
    if(lambda.min==0){
      lambda <- exp(seq(log(0.001*lambda.max), log(lambda.max), length.out = nlambda))
    }else{
      lambda <- exp(seq(log(lambda.min*lambda.max), log(lambda.max), length.out = nlambda))
    }
  }
  beta_lambda <- matrix(0, nrow = p, ncol = length(lambda))
  beta <- rep(0, p) -> beta_new
  iter <- rep(0, length(lambda))
  z <- rep(0, p)
  if(penalty!="lasso"){
    for (i in 1:length(lambda)) {
      r <- y - X%*%beta
      t <- 0
      repeat{
        for (j in 1:p) {
          z[j] <- X[, j]%*%r/n + beta[j]
          beta_new[j] <- get(penalty)(z = z[j], lambda = lambda[i], gamma = gamma)
          r <- r - (beta_new[j]-beta[j])*X[, j]
        }
        t <- t + 1
        if(norm(beta_new-beta,"2")<epsilon){break}
        beta <- beta_new
      }
      beta_lambda[, i] <- beta
      iter[i] <- t
    }
  }else{
    for (i in 1:length(lambda)) {
      r <- y - X%*%beta
      t <- 0
      repeat{
        for (j in 1:p) {
          z[j] <- X[, j]%*%r/n + beta[j]
          beta_new[j] <- S(z = z[j], lambda = lambda[i])
          r <- r - (beta_new[j]-beta[j])*X[, j]
        }
        t <- t + 1
        if(norm(beta_new-beta,"2")<epsilon){break}
        beta <- beta_new
      }
      beta_lambda[, i] <- beta
      iter[i] <- t
    }
  }
  
  if(plot){
    if(penalty!="lasso"){
      for (i in length(lambda):1) {
        if(i==1){
          U <- which(beta_lambda[, i]!=0)
        }else{
          U <- which(beta_lambda[, i]!=0|beta_lambda[, i-1]!=0)
        }
        c_star <- eigen(t(X[, U])%*%X[, U]/n)$values[length(U)]
        if(gamma<=1/c_star){break}
      }
      lambda_star <- lambda[i]
      
      plot(lambda, beta_lambda[1, ], "l", col = rgb(runif(1), runif(1), runif(1)), 
        ylim = c(min(beta_lambda), max(beta_lambda)), xlab = expression(lambda), 
        ylab = expression(hat(beta)), xlim = c(lambda[length(lambda)], lambda[1]))
      apply(matrix(2:p), 1, function(x){
        lines(lambda, beta_lambda[x, ], col = rgb(runif(1), runif(1), runif(1)))
      })
      abline(h = 0)
      polygon(x = c(lambda_star, lambda_star, lambda[1], lambda[1]), 
        y = c(min(beta_lambda), max(beta_lambda), max(beta_lambda), min(beta_lambda)), 
        col = gray(0.5, 0.5), border = NA)
    }else{
      plot(lambda, beta_lambda[1, ], "l", col = rgb(runif(1), runif(1), runif(1)), 
        ylim = c(min(beta_lambda), max(beta_lambda)), xlab = expression(lambda), 
        ylab = expression(hat(beta)), xlim = c(lambda[length(lambda)], lambda[1]))
      apply(matrix(2:p), 1, function(x){
        lines(lambda, beta_lambda[x, ], col = rgb(runif(1), runif(1), runif(1)))
      })
      abline(h = 0)
    }
    
  }

  
  return(list(beta = beta_lambda, lambda = lambda, iter = iter))
}
