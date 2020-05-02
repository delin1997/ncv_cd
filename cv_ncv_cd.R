cv_ncv_cd <- function(X, y, family=c("gaussian", "binomial"), lambda.min, lambda, nlambda=100, penalty=c("MCP", "SCAD", "lasso"), gamma=switch(penalty, SCAD=3.7, MCP=3), nfolds=10){
  n <- nrow(X)
  index <- sample(1:nfolds, replace = T, size = n)
  correct <- rep(0, )
}